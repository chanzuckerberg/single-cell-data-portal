import logging
from functools import wraps
from urllib.parse import unquote

from flask import (
    Blueprint,
    current_app,
    redirect,
    request,
    send_from_directory,
)
from flask_restful import Api, Resource

import server.common.rest as common_rest
from server.app.api.util import get_data_adaptor, get_dataset_artifact_s3_uri
from server.common.errors import (
    DatasetAccessError,
    DatasetMetadataError,
    DatasetNotFoundError,
    TombstoneError,
)
from server.common.utils.http_cache import ONE_YEAR, cache_control


def rest_get_s3uri_data_adaptor(func):
    @wraps(func)
    def wrapped_function(self, s3_uri=None):
        try:
            s3_uri = unquote(s3_uri) if s3_uri else s3_uri
            data_adaptor = get_data_adaptor(s3_uri)
            return func(self, data_adaptor)
        except (DatasetAccessError, DatasetNotFoundError, DatasetMetadataError) as e:
            return common_rest.abort_and_log(
                e.status_code, f"Invalid s3_uri {s3_uri}: {e.message}", loglevel=logging.INFO, include_exc_info=True
            )

    return wrapped_function


class DatasetResource(Resource):
    """Base class for all Resources that act on datasets."""

    def __init__(self, url_dataroot):
        super().__init__()
        self.url_dataroot = url_dataroot


class S3URIAPI(DatasetResource):
    @cache_control(public=True, no_store=True, max_age=0)
    def get(self, dataset):
        app_config = current_app.app_config
        return common_rest.s3_uri_get(app_config, self.url_dataroot, dataset)


class S3URIResource(Resource):
    """Base class for all Resources that act on S3 URIs."""

    def __init__(self, s3_uri_root):
        super().__init__()
        self.s3_uri_root = s3_uri_root


class SchemaAPI(S3URIResource):
    # TODO @mdunitz separate dataset schema and user schema
    @cache_control(immutable=True, max_age=ONE_YEAR)
    @rest_get_s3uri_data_adaptor
    def get(self, data_adaptor):
        return common_rest.schema_get(data_adaptor)


class ConfigAPI(S3URIResource):
    @cache_control(immutable=True, max_age=ONE_YEAR)
    @rest_get_s3uri_data_adaptor
    def get(self, data_adaptor):
        return common_rest.config_get(current_app.app_config, data_adaptor)


class AnnotationsObsAPI(S3URIResource):
    @cache_control(immutable=True, max_age=ONE_YEAR)
    @rest_get_s3uri_data_adaptor
    def get(self, data_adaptor):
        return common_rest.annotations_obs_get(request, data_adaptor)


class AnnotationsVarAPI(S3URIResource):
    @cache_control(immutable=True, max_age=ONE_YEAR)
    @rest_get_s3uri_data_adaptor
    def get(self, data_adaptor):
        return common_rest.annotations_var_get(request, data_adaptor)


class DataVarAPI(S3URIResource):
    @cache_control(no_store=True)
    @rest_get_s3uri_data_adaptor
    def put(self, data_adaptor):
        return common_rest.data_var_put(request, data_adaptor)

    @cache_control(immutable=True, max_age=ONE_YEAR)
    @rest_get_s3uri_data_adaptor
    def get(self, data_adaptor):
        return common_rest.data_var_get(request, data_adaptor)


class ColorsAPI(S3URIResource):
    @cache_control(immutable=True, max_age=ONE_YEAR)
    @rest_get_s3uri_data_adaptor
    def get(self, data_adaptor):
        return common_rest.colors_get(data_adaptor)


class DiffExpObsAPI(DatasetResource):
    @cache_control(no_store=True)
    @rest_get_s3uri_data_adaptor
    def post(self, data_adaptor):
        return common_rest.diffexp_obs_post(request, data_adaptor)


class DiffExpObs2API(DatasetResource):
    @cache_control(no_store=True)
    @rest_get_s3uri_data_adaptor
    def post(self, data_adaptor):
        return common_rest.diffex_binary_post(request, data_adaptor)


class LayoutObsAPI(S3URIResource):
    @cache_control(immutable=True, max_age=ONE_YEAR)
    @rest_get_s3uri_data_adaptor
    def get(self, data_adaptor):
        return common_rest.layout_obs_get(request, data_adaptor)


class SummarizeVarAPI(S3URIResource):
    @rest_get_s3uri_data_adaptor
    @cache_control(immutable=True, max_age=ONE_YEAR)
    def get(self, data_adaptor):
        return common_rest.summarize_var_get(request, data_adaptor)

    @rest_get_s3uri_data_adaptor
    @cache_control(no_store=True)
    def post(self, data_adaptor):
        return common_rest.summarize_var_post(request, data_adaptor)


class GeneInfoAPI(S3URIResource):
    @rest_get_s3uri_data_adaptor
    def get(self, data_adaptor):
        return common_rest.gene_info_get(request)


class VersionAPI(Resource):
    def get(self):
        return common_rest.get_deployed_version(request)


def rest_get_dataset_explorer_location_data_adaptor(func):
    @wraps(func)
    def wrapped_function(self, dataset=None):
        try:
            s3_uri = get_dataset_artifact_s3_uri(self.url_dataroot, dataset)
            data_adaptor = get_data_adaptor(s3_uri)
            # HACK: Used *only* to pass the dataset_explorer_location to DatasetMeta.get_dataset_and_collection_
            # metadata()
            data_adaptor.dataset_id = dataset
            return func(self, data_adaptor)
        except (DatasetAccessError, DatasetNotFoundError, DatasetMetadataError) as e:
            return common_rest.abort_and_log(
                e.status_code, f"Invalid s3_uri {dataset}: {e.message}", loglevel=logging.INFO, include_exc_info=True
            )
        except TombstoneError as e:
            parent_collection_url = (
                f"{current_app.app_config.server__app__web_base_url}/collections/{e.collection_id}"  # noqa E501
            )
            return redirect(f"{parent_collection_url}?tombstoned_dataset_id={e.dataset_id}")

    return wrapped_function


class DatasetMetadataAPI(DatasetResource):
    @cache_control(public=True, no_store=True, max_age=0)
    @rest_get_dataset_explorer_location_data_adaptor
    def get(self, data_adaptor):
        return common_rest.dataset_metadata_get(current_app.app_config, self.url_dataroot, data_adaptor.dataset_id)


def get_api_dataroot_resources(bp_dataroot, url_dataroot=None):
    """Add resources that refer to a dataset"""
    api = Api(bp_dataroot)

    def add_resource(resource, url):
        """convenience function to make the outer function less verbose"""
        api.add_resource(resource, url, resource_class_args=(url_dataroot,))

    # Initialization routes
    add_resource(S3URIAPI, "/s3_uri")
    add_resource(DatasetMetadataAPI, "/dataset-metadata")
    return api


def get_api_s3uri_resources(bp_dataroot, s3uri_path):
    """Add resources that refer to a S3 URIs"""
    api = Api(bp_dataroot)

    def add_resource(resource, url):
        """convenience function to make the outer function less verbose"""
        api.add_resource(resource, url, resource_class_args=(s3uri_path,))

    # Initialization routes
    add_resource(SchemaAPI, "/schema")
    add_resource(ConfigAPI, "/config")
    # Data routes
    add_resource(AnnotationsObsAPI, "/annotations/obs")
    add_resource(AnnotationsVarAPI, "/annotations/var")
    add_resource(DataVarAPI, "/data/var")
    add_resource(GeneInfoAPI, "/geneinfo")
    add_resource(SummarizeVarAPI, "/summarize/var")
    # Display routes
    add_resource(ColorsAPI, "/colors")
    # Computation routes
    add_resource(DiffExpObsAPI, "/diffexp/obs")
    add_resource(DiffExpObs2API, "/diffexp/obs2")
    add_resource(LayoutObsAPI, "/layout/obs")
    return api


def register_api_v3(app, app_config, api_url_prefix):
    api_version = "/api/v0.3"

    s3uri_api_path = "s3_uri"
    bp_s3uri = Blueprint(
        f"api_dataset_{s3uri_api_path}_{api_version.replace('.',',')}",
        __name__,
        url_prefix=(f"{api_url_prefix}/{s3uri_api_path}/<s3_uri>" + api_version).replace("//", "/"),
    )
    s3uri_resources = get_api_s3uri_resources(bp_s3uri, s3uri_api_path)
    app.register_blueprint(s3uri_resources.blueprint)

    Api(app).add_resource(VersionAPI, "/deployed_version")

    # NOTE:  These routes only allow the dataset to be in the directory
    # of the dataroot, and not a subdirectory.  We may want to change
    # the route format at some point
    for dataroot_dict in app_config.server__multi_dataset__dataroots.values():
        url_dataroot = dataroot_dict["base_url"]
        bp_dataroot = Blueprint(
            name=f"api_dataset_{url_dataroot}_{api_version.replace('.',',')}",
            import_name=__name__,
            url_prefix=(f"{api_url_prefix}/{url_dataroot}/<string:dataset>" + api_version).replace("//", "/"),
        )
        dataroot_resources = get_api_dataroot_resources(bp_dataroot, url_dataroot)
        app.register_blueprint(dataroot_resources.blueprint)
        app.add_url_rule(
            f"/{url_dataroot}/<string:dataset>/static/<path:filename>",
            f"static_assets_{url_dataroot}",
            view_func=lambda dataset, filename: send_from_directory("../common/web/static", filename),
            methods=["GET"],
        )
