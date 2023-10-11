import copy
import hashlib
import logging
import os
import struct
import sys
import zlib
from http import HTTPStatus

import requests
from flask import abort, current_app, jsonify, make_response, redirect
from werkzeug.urls import url_unquote

from server.app.api.util import get_dataset_artifact_s3_uri
from server.common.config.client_config import get_client_config
from server.common.constants import Axis, DiffExpMode, JSON_NaN_to_num_warning_msg
from server.common.diffexpdu import DiffExArguments
from server.common.errors import (
    ColorFormatException,
    DatasetAccessError,
    DisabledFeatureError,
    ExceedsLimitError,
    FilterError,
    InvalidCxgDatasetError,
    JSONEncodingValueError,
    TombstoneError,
    UnsupportedSummaryMethod,
)
from server.dataset import dataset_metadata


def abort_and_log(code, logmsg, loglevel=logging.DEBUG, include_exc_info=False):
    """
    Log the message, then abort with HTTP code. If include_exc_info is true,
    also include current exception via sys.exc_info().
    """
    exc_info = sys.exc_info() if include_exc_info else False
    current_app.logger.log(loglevel, logmsg, exc_info=exc_info)
    # Do NOT send log message to HTTP response.
    return abort(code)


def _query_parameter_to_filter(args):
    """
    Convert an annotation value filter, if present in the query args,
    into the standard dict filter format used by internal code.

    Query param filters look like:  <axis>:name=value, where value
    may be one of:
        - a range, min,max, where either may be an open range by using an asterisc, eg, 10,*
        - a value
    Eg,
        ...?tissue=lung&obs:tissue=heart&obs:num_reads=1000,*
    """
    filters = {
        "obs": {},
        "var": {},
    }

    # args has already been url-unquoted once.  We assume double escaping
    # on name and value.
    try:
        for key, value in args.items(multi=True):
            axis, name = key.split(":")
            if axis not in ("obs", "var"):
                raise FilterError("unknown filter axis")
            name = url_unquote(name)
            current = filters[axis].setdefault(name, {"name": name})

            val_split = value.split(",")
            if len(val_split) == 1:
                if "min" in current or "max" in current:
                    raise FilterError("do not mix range and value filters")
                value = url_unquote(value)
                values = current.setdefault("values", [])
                values.append(value)

            elif len(val_split) == 2:
                if len(current) > 1:
                    raise FilterError("duplicate range specification")
                min = url_unquote(val_split[0])
                max = url_unquote(val_split[1])
                if min != "*":
                    current["min"] = float(min)
                if max != "*":
                    current["max"] = float(max)
                if len(current) < 2:
                    raise FilterError("must specify at least min or max in range filter")

            else:
                raise FilterError("badly formated filter value")

    except ValueError as e:
        raise FilterError(str(e)) from None

    result = {}
    for axis in ("obs", "var"):
        axis_filter = filters[axis]
        if len(axis_filter) > 0:
            result[axis] = {"annotation_value": list(axis_filter.values())}

    return result


def schema_get_helper(data_adaptor):
    """helper function to gather the schema from the data source and annotations"""
    schema = data_adaptor.get_schema()
    schema = copy.deepcopy(schema)

    return schema


def schema_get(data_adaptor):
    schema = schema_get_helper(data_adaptor)
    return make_response(jsonify({"schema": schema}), HTTPStatus.OK)


def dataset_metadata_get(app_config, url_dataroot, dataset_id):
    metadata = dataset_metadata.get_dataset_and_collection_metadata(url_dataroot, dataset_id, app_config)
    if metadata is not None:
        return make_response(jsonify({"metadata": metadata}), HTTPStatus.OK)
    else:
        return abort(HTTPStatus.NOT_FOUND)


def s3_uri_get(app_config, url_dataroot_id, dataset_id):
    try:
        dataset_artifact_s3_uri = get_dataset_artifact_s3_uri(url_dataroot_id, dataset_id)
    except TombstoneError as e:
        parent_collection_url = (
            f"{current_app.app_config.server__app__web_base_url}/collections/{e.collection_id}"  # noqa E501
        )
        return redirect(f"{parent_collection_url}?tombstoned_dataset_id={e.dataset_id}")
    else:
        return make_response(jsonify(dataset_artifact_s3_uri), HTTPStatus.OK)


def config_get(app_config, data_adaptor):
    config = get_client_config(app_config, data_adaptor, current_app)
    return make_response(jsonify(config), HTTPStatus.OK)


def annotations_obs_get(request, data_adaptor):
    fields = request.args.getlist("annotation-name", None)
    nBins = request.args.get("nbins", None)
    if nBins is not None:
        nBins = int(nBins)

    num_columns_requested = len(data_adaptor.get_obs_keys()) if len(fields) == 0 else len(fields)
    if data_adaptor.app_config.exceeds_limit("column_request_max", num_columns_requested):
        return abort(HTTPStatus.BAD_REQUEST)
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    try:
        fbs = data_adaptor.annotation_to_fbs_matrix(Axis.OBS, fields, num_bins=nBins)
        return make_response(fbs, HTTPStatus.OK, {"Content-Type": "application/octet-stream"})
    except KeyError as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def inflate(data):
    return zlib.decompress(data)


def annotations_var_get(request, data_adaptor):
    fields = request.args.getlist("annotation-name", None)
    nBins = request.args.get("nbins", None)
    if nBins is not None:
        nBins = int(nBins)

    num_columns_requested = len(data_adaptor.get_var_keys()) if len(fields) == 0 else len(fields)
    if data_adaptor.app_config.exceeds_limit("column_request_max", num_columns_requested):
        return abort(HTTPStatus.BAD_REQUEST)
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    try:
        return make_response(
            data_adaptor.annotation_to_fbs_matrix(Axis.VAR, fields, num_bins=nBins),
            HTTPStatus.OK,
            {"Content-Type": "application/octet-stream"},
        )
    except KeyError as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def data_var_put(request, data_adaptor):
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    filter_json = request.get_json()
    filter = filter_json["filter"] if filter_json else None
    nBins = filter_json.get("nbins", None) if filter_json else None
    if nBins is not None:
        nBins = int(nBins)

    try:
        return make_response(
            data_adaptor.data_frame_to_fbs_matrix(filter, axis=Axis.VAR, num_bins=nBins),
            HTTPStatus.OK,
            {"Content-Type": "application/octet-stream"},
        )
    except (FilterError, ValueError, ExceedsLimitError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def data_var_get(request, data_adaptor):
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    try:
        nBins = request.args.get("nbins", None)
        if nBins is not None:
            nBins = int(nBins)
        args_filter_only = request.args.copy()
        args_filter_only.poplist("nbins")
        filter = _query_parameter_to_filter(args_filter_only)
        return make_response(
            data_adaptor.data_frame_to_fbs_matrix(filter, axis=Axis.VAR, num_bins=nBins),
            HTTPStatus.OK,
            {"Content-Type": "application/octet-stream"},
        )
    except (FilterError, ValueError, ExceedsLimitError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def colors_get(data_adaptor):
    if not data_adaptor.app_config.default_dataset__presentation__custom_colors:
        return make_response(jsonify({}), HTTPStatus.OK)
    try:
        return make_response(jsonify(data_adaptor.get_colors()), HTTPStatus.OK)
    except ColorFormatException as e:
        return abort_and_log(HTTPStatus.NOT_FOUND, str(e), include_exc_info=True)


def gene_info_get(request):
    """
    Request information about a gene from the data portal gene_info api
    """
    api_base_url = current_app.app_config.server__gene_info__api_base
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    try:
        response = requests.get(
            url=f"{api_base_url}/gene_info?geneID={request.args['geneID']}&gene={request.args['gene']}", headers=headers
        )
        if response.status_code == 200:
            return make_response(response.content, HTTPStatus.OK, {"Content-Type": "application/json"})
        else:
            # in the event of a failed search, return empty response
            return None
    except Exception as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)


def get_deployed_version(request):
    """
    Returns the deployed version
    """
    version_info = {"Explorer": os.getenv("COMMIT_SHA")}
    return make_response(jsonify(version_info), 200)


def diffexp_obs_post(request, data_adaptor):
    if not data_adaptor.app_config.default_dataset__diffexp__enable:
        return abort(HTTPStatus.NOT_IMPLEMENTED)

    args = request.get_json()
    try:
        # TODO: implement varfilter mode
        mode = DiffExpMode(args["mode"])
        if mode == DiffExpMode.VAR_FILTER or "varFilter" in args:
            return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, "varFilter not enabled")

        set1_filter = args.get("set1", {"filter": {}})["filter"]
        set2_filter = args.get("set2", {"filter": {}})["filter"]
        count = data_adaptor.app_config.default_dataset__diffexp__count

        if set1_filter is None or set2_filter is None or count is None:
            return abort_and_log(HTTPStatus.BAD_REQUEST, "missing required parameter")
        if Axis.VAR in set1_filter or Axis.VAR in set2_filter:
            return abort_and_log(HTTPStatus.BAD_REQUEST, "var axis filter not enabled")

    except (KeyError, TypeError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)

    try:
        diffexp = data_adaptor.diffexp_topN(set1_filter, set2_filter, count)
        return make_response(diffexp, HTTPStatus.OK, {"Content-Type": "application/json"})
    except (ValueError, DisabledFeatureError, FilterError, ExceedsLimitError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)
    except JSONEncodingValueError:
        # JSON encoding failure, usually due to bad data. Just let it ripple up
        # to default exception handler.
        current_app.logger.warning(JSON_NaN_to_num_warning_msg)
        raise


def diffex_binary_post(request, data_adaptor):
    MAX_CONTENT_LENGTH = 100 * 1024**2  # perhaps a config variable would be suitable?
    if not data_adaptor.app_config.default_dataset__diffexp__enable:
        return abort(HTTPStatus.NOT_IMPLEMENTED)
    if not request.content_type or "application/octet-stream" not in request.content_type:
        return abort(HTTPStatus.UNSUPPORTED_MEDIA_TYPE)
    if not request.content_length:
        return abort(HTTPStatus.LENGTH_REQUIRED)
    if request.content_length > MAX_CONTENT_LENGTH:
        return abort(HTTPStatus.REQUEST_ENTITY_TOO_LARGE)

    try:
        buf = request.get_data()
        diffex_args = DiffExArguments.unpack_from(buf)

        if diffex_args.mode != DiffExArguments.DiffExMode.TopN:
            return abort_and_log(HTTPStatus.NOT_IMPLEMENTED, "Unsupported diffex mode")

        result = data_adaptor.diffexp_topN_from_list(diffex_args.set1, diffex_args.set2, diffex_args.params.N)
        return make_response(result, HTTPStatus.OK, {"Content-Type": "application/json"})

    except (KeyError, TypeError, AssertionError, struct.error) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)
    except JSONEncodingValueError:
        # JSON encoding failure, usually due to bad data. Just let it ripple up
        # to default exception handler.
        current_app.logger.warning(JSON_NaN_to_num_warning_msg)
        raise


def layout_obs_get(request, data_adaptor):
    fields = request.args.getlist("layout-name", None)
    nBins = request.args.get("nbins", None)
    if nBins is not None:
        nBins = int(nBins)

    num_columns_requested = len(data_adaptor.get_embedding_names()) if len(fields) == 0 else len(fields)
    if data_adaptor.app_config.exceeds_limit("column_request_max", num_columns_requested):
        return abort(HTTPStatus.BAD_REQUEST)

    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    try:
        return make_response(
            data_adaptor.layout_to_fbs_matrix(fields, num_bins=nBins),
            HTTPStatus.OK,
            {"Content-Type": "application/octet-stream"},
        )
    except (KeyError, DatasetAccessError) as e:
        return abort_and_log(HTTPStatus.BAD_REQUEST, str(e), include_exc_info=True)
    except InvalidCxgDatasetError:
        return abort_and_log(
            HTTPStatus.NOT_IMPLEMENTED,
            f"No embedding available {request.path}",
            loglevel=logging.ERROR,
            include_exc_info=True,
        )


def summarize_var_helper(request, data_adaptor, key, raw_query):
    preferred_mimetype = request.accept_mimetypes.best_match(["application/octet-stream"])
    if preferred_mimetype != "application/octet-stream":
        return abort(HTTPStatus.NOT_ACCEPTABLE)

    summary_method = request.values.get("method", default="mean")
    nBins = request.values.get("nbins", default=None)
    if nBins is not None:
        nBins = int(nBins)

    def summarizeQueryHash(raw_query):
        """generate a cache key (hash) from the raw query string"""
        return hashlib.sha1(raw_query).hexdigest()

    query_hash = summarizeQueryHash(raw_query)
    if key and query_hash != key:
        return abort(HTTPStatus.BAD_REQUEST, description="query key did not match")

    args_filter_only = request.values.copy()
    args_filter_only.poplist("method")
    args_filter_only.poplist("key")
    args_filter_only.poplist("nbins")

    try:
        filter = _query_parameter_to_filter(args_filter_only)
        return make_response(
            data_adaptor.summarize_var(summary_method, filter, query_hash, num_bins=nBins),
            HTTPStatus.OK,
            {"Content-Type": "application/octet-stream"},
        )
    except ValueError as e:
        return abort(HTTPStatus.NOT_FOUND, description=str(e))
    except (UnsupportedSummaryMethod, FilterError) as e:
        return abort(HTTPStatus.BAD_REQUEST, description=str(e))


def summarize_var_get(request, data_adaptor):
    return summarize_var_helper(request, data_adaptor, None, request.query_string)


def summarize_var_post(request, data_adaptor):
    if not request.content_type or "application/x-www-form-urlencoded" not in request.content_type:
        return abort(HTTPStatus.UNSUPPORTED_MEDIA_TYPE)
    if request.content_length > 1_000_000:  # just a sanity check to avoid memory exhaustion
        return abort(HTTPStatus.BAD_REQUEST)

    key = request.args.get("key", default=None)
    return summarize_var_helper(request, data_adaptor, key, request.get_data())
