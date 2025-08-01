# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
import re
import pyspark.pandas as ps
from pyspark.pandas.internal import InternalFrame
from pyspark.pandas.series import first_series
from pyspark.pandas.utils import scol_for
from pyspark.sql.functions import expr, col, lit
from pyspark.sql.types import StructType, StructField, StringType, IntegerType

from sedona.geopandas import GeoDataFrame, GeoSeries

# Pre-compiled regex pattern for suffix validation
SUFFIX_PATTERN = re.compile(r"^[a-zA-Z_][a-zA-Z0-9_]*$")


def _frame_join(
    left_df,
    right_df,
    how="inner",
    predicate="intersects",
    lsuffix="left",
    rsuffix="right",
    distance=None,
    on_attribute=None,
):
    """Join the GeoDataFrames at the DataFrame level.

    Parameters
    ----------
    left_df : GeoDataFrame or GeoSeries
        Left dataset to join
    right_df : GeoDataFrame or GeoSeries
        Right dataset to join
    how : str, default 'inner'
        Join type: 'inner', 'left', 'right'
    predicate : str, default 'intersects'
        Spatial predicate to use
    lsuffix : str, default 'left'
        Suffix for left overlapping columns
    rsuffix : str, default 'right'
        Suffix for right overlapping columns
    distance : float, optional
        Distance parameter for dwithin predicate
    on_attribute : list, optional
        Additional columns to join on

    Returns
    -------
    GeoDataFrame or GeoSeries
        Joined result
    """
    # Predicate mapping
    predicate_map = {
        "intersects": "ST_Intersects",
        "contains": "ST_Contains",
        "within": "ST_Within",
        "touches": "ST_Touches",
        "crosses": "ST_Crosses",
        "overlaps": "ST_Overlaps",
        "dwithin": "ST_DWithin",
    }

    if predicate not in predicate_map:
        raise ValueError(
            f"Predicate '{predicate}' not supported. Available: {list(predicate_map.keys())}"
        )

    spatial_func = predicate_map[predicate]

    # Get the internal Spark DataFrames
    left_sdf = left_df._internal.spark_frame
    right_sdf = right_df._internal.spark_frame

    # Handle geometry columns - check if they exist and get proper column names
    left_geom_col = None
    right_geom_col = None

    # Find geometry columns in left dataframe
    for field in left_sdf.schema.fields:
        if field.dataType.typeName() in ("geometrytype", "binary"):
            left_geom_col = field.name
            break

    # Find geometry columns in right dataframe
    for field in right_sdf.schema.fields:
        if field.dataType.typeName() in ("geometrytype", "binary"):
            right_geom_col = field.name
            break

    if left_geom_col is None or right_geom_col is None:
        raise ValueError("Both datasets must have geometry columns")

    # Prepare geometry expressions for spatial join
    if left_sdf.schema[left_geom_col].dataType.typeName() == "binary":
        left_geom_expr = f"ST_GeomFromWKB(`{left_geom_col}`) as l_geometry"
    else:
        left_geom_expr = f"`{left_geom_col}` as l_geometry"

    if right_sdf.schema[right_geom_col].dataType.typeName() == "binary":
        right_geom_expr = f"ST_GeomFromWKB(`{right_geom_col}`) as r_geometry"
    else:
        right_geom_expr = f"`{right_geom_col}` as r_geometry"

    # Select all columns with geometry
    left_cols = [left_geom_expr] + [
        f"`{field.name}` as l_{field.name}"
        for field in left_sdf.schema.fields
        if field.name != left_geom_col and not field.name.startswith("__")
    ]
    right_cols = [right_geom_expr] + [
        f"`{field.name}` as r_{field.name}"
        for field in right_sdf.schema.fields
        if field.name != right_geom_col and not field.name.startswith("__")
    ]

    left_geo_df = left_sdf.selectExpr(*left_cols)
    right_geo_df = right_sdf.selectExpr(*right_cols)

    # Build spatial join condition
    if predicate == "dwithin":
        if distance is None:
            raise ValueError("Distance parameter is required for 'dwithin' predicate")
        spatial_condition = f"{spatial_func}(l_geometry, r_geometry, {distance})"
    else:
        spatial_condition = f"{spatial_func}(l_geometry, r_geometry)"

    # Add attribute-based join condition if specified
    join_condition = spatial_condition
    if on_attribute:
        for attr in on_attribute:
            join_condition += f" AND l_{attr} = r_{attr}"

    # Perform spatial join based on join type
    if how == "inner":
        spatial_join_df = left_geo_df.alias("l").join(
            right_geo_df.alias("r"), expr(join_condition)
        )
    elif how == "left":
        spatial_join_df = left_geo_df.alias("l").join(
            right_geo_df.alias("r"), expr(join_condition), "left"
        )
    elif how == "right":
        spatial_join_df = left_geo_df.alias("l").join(
            right_geo_df.alias("r"), expr(join_condition), "right"
        )
    else:
        raise ValueError(f"Join type '{how}' not supported")

    # Handle column naming with suffixes
    final_columns = []

    # Add geometry column (always from left for geopandas compatibility)
    # Currently, Sedona stores geometries in EWKB format
    final_columns.append("ST_AsEWKB(l_geometry) as geometry")

    # Add other columns with suffix handling
    left_data_cols = [col for col in left_geo_df.columns if col != "l_geometry"]
    right_data_cols = [col for col in right_geo_df.columns if col != "r_geometry"]

    for col_name in left_data_cols:
        base_name = col_name[2:]  # Remove "l_" prefix
        right_col = f"r_{base_name}"

        if right_col in right_data_cols:
            # Column exists in both - apply suffixes
            final_columns.append(f"{col_name} as {base_name}_{lsuffix}")
        else:
            # Column only in left
            final_columns.append(f"{col_name} as {base_name}")

    for col_name in right_data_cols:
        base_name = col_name[2:]  # Remove "r_" prefix
        left_col = f"l_{base_name}"

        if left_col in left_data_cols:
            # Column exists in both - apply suffixes
            final_columns.append(f"{col_name} as {base_name}_{rsuffix}")
        else:
            # Column only in right
            final_columns.append(f"{col_name} as {base_name}")

    # Select final columns
    result_df = spatial_join_df.selectExpr(*final_columns)

    # Return appropriate type based on input
    if isinstance(left_df, GeoSeries) and isinstance(right_df, GeoSeries):
        # Return GeoSeries for GeoSeries inputs
        internal = InternalFrame(
            spark_frame=result_df,
            index_spark_columns=None,
            column_labels=[left_df._col_label],
            data_spark_columns=[scol_for(result_df, "geometry")],
            data_fields=[left_df._internal.data_fields[0]],
            column_label_names=left_df._internal.column_label_names,
        )
        return _to_geo_series(first_series(ps.DataFrame(internal)))
    else:
        # Return GeoDataFrame for GeoDataFrame inputs
        return GeoDataFrame(result_df)


def sjoin(
    left_df,
    right_df,
    how="inner",
    predicate="intersects",
    lsuffix="left",
    rsuffix="right",
    distance=None,
    on_attribute=None,
    **kwargs,
):
    """Spatial join of two GeoDataFrames.

    Parameters
    ----------
    left_df, right_df : GeoDataFrames
    how : string, default 'inner'
        The type of join:

        * 'left': use keys from left_df; retain only left_df geometry column
        * 'right': use keys from right_df; retain only right_df geometry column
        * 'inner': use intersection of keys from both dfs; retain only
          left_df geometry column
    predicate : string, default 'intersects'
        Binary predicate. Valid values are determined by the spatial index used.
        You can check the valid values in left_df or right_df as
        ``left_df.sindex.valid_query_predicates`` or
        ``right_df.sindex.valid_query_predicates``
        Replaces deprecated ``op`` parameter.
    lsuffix : string, default 'left'
        Suffix to apply to overlapping column names (left GeoDataFrame).
    rsuffix : string, default 'right'
        Suffix to apply to overlapping column names (right GeoDataFrame).
    distance : number or array_like, optional
        Distance(s) around each input geometry within which to query the tree
        for the 'dwithin' predicate. If array_like, must be
        one-dimesional with length equal to length of left GeoDataFrame.
        Required if ``predicate='dwithin'``.
    on_attribute : string, list or tuple
        Column name(s) to join on as an additional join restriction on top
        of the spatial predicate. These must be found in both DataFrames.
        If set, observations are joined only if the predicate applies
        and values in specified columns match.

    Examples
    --------
    >>> groceries_w_communities = geopandas.sjoin(groceries, chicago)
    >>> groceries_w_communities.head()  # doctest: +SKIP
       OBJECTID       community                           geometry
    0        16          UPTOWN  MULTIPOINT ((-87.65661 41.97321))
    1        18     MORGAN PARK  MULTIPOINT ((-87.68136 41.69713))
    2        22  NEAR WEST SIDE  MULTIPOINT ((-87.63918 41.86847))
    3        23  NEAR WEST SIDE  MULTIPOINT ((-87.65495 41.87783))
    4        27         CHATHAM  MULTIPOINT ((-87.62715 41.73623))
    [5 rows x 95 columns]

    Notes
    -----
    Every operation in GeoPandas is planar, i.e. the potential third
    dimension is not taken into account.
    """
    if kwargs:
        first = next(iter(kwargs.keys()))
        raise TypeError(f"sjoin() got an unexpected keyword argument '{first}'")

    on_attribute = _maybe_make_list(on_attribute)

    _basic_checks(left_df, right_df, how, lsuffix, rsuffix, on_attribute=on_attribute)

    joined = _frame_join(
        left_df,
        right_df,
        how=how,
        predicate=predicate,
        lsuffix=lsuffix,
        rsuffix=rsuffix,
        distance=distance,
        on_attribute=on_attribute,
    )

    return joined


def _maybe_make_list(obj):
    if isinstance(obj, tuple):
        return list(obj)
    if obj is not None and not isinstance(obj, list):
        return [obj]
    return obj


def _basic_checks(left_df, right_df, how, lsuffix, rsuffix, on_attribute=None):
    """Checks the validity of join input parameters.

    `how` must be one of the valid options.
    `'index_'` concatenated with `lsuffix` or `rsuffix` must not already
    exist as columns in the left or right data frames.

    Parameters
    ------------
    left_df : GeoDataFrame or GeoSeries
    right_df : GeoDataFrame or GeoSeries
    how : str, one of 'left', 'right', 'inner'
        join type
    lsuffix : str
        left index suffix
    rsuffix : str
        right index suffix
    on_attribute : list, default None
        list of column names to merge on along with geometry
    """
    if not isinstance(left_df, (GeoSeries, GeoDataFrame)):
        raise ValueError(
            f"'left_df' should be GeoSeries or GeoDataFrame, got {type(left_df)}"
        )

    if not isinstance(right_df, (GeoSeries, GeoDataFrame)):
        raise ValueError(
            f"'right_df' should be GeoSeries or GeoDataFrame, got {type(right_df)}"
        )

    allowed_hows = ["inner", "left", "right"]
    if how not in allowed_hows:
        raise ValueError(f'`how` was "{how}" but is expected to be in {allowed_hows}')

    # Check if on_attribute columns exist in both datasets
    if on_attribute:
        for attr in on_attribute:
            if hasattr(left_df, "columns") and attr not in left_df.columns:
                raise ValueError(f"Column '{attr}' not found in left dataset")
            if hasattr(right_df, "columns") and attr not in right_df.columns:
                raise ValueError(f"Column '{attr}' not found in right dataset")

    # Check for reserved column names that would conflict
    if lsuffix == rsuffix:
        raise ValueError("lsuffix and rsuffix cannot be the same")

    # Validate suffix format (should not contain special characters that would break SQL)
    if not SUFFIX_PATTERN.match(lsuffix):
        raise ValueError(f"lsuffix '{lsuffix}' contains invalid characters")
    if not SUFFIX_PATTERN.match(rsuffix):
        raise ValueError(f"rsuffix '{rsuffix}' contains invalid characters")


def _to_geo_series(df: ps.Series) -> GeoSeries:
    """
    Get the first Series from the DataFrame.

    Parameters:
    - df: The input DataFrame.

    Returns:
    - GeoSeries: The first Series from the DataFrame.
    """
    return GeoSeries(data=df)
