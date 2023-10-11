"""
This module contains the configuration model used to validate the cellxgene configuration.

Note: `root_validator`s and class variables are executed in the order they are defined in the class definition.
"""

from __future__ import annotations

import os
import sys
import warnings
from enum import Enum
from typing import Dict, List, Optional, Union
from urllib.parse import quote_plus

from pydantic import BaseModel, Extra, Field, root_validator, validator

from server.common.utils.data_locator import discover_s3_region_name
from server.common.utils.utils import custom_format_warning, find_available_port, is_port_available


class CspDirectives(BaseModel):
    img_src: Union[str, List[str]] = Field(default_factory=list, alias="img-src")
    script_src: Union[str, List[str]] = Field(default_factory=list, alias="script-src")
    connect_src: Union[str, List[str]] = Field(default_factory=list, alias="connect-src")

    @root_validator(skip_on_failure=True)
    def string_to_list(cls, values):
        for key, value in values.items():
            if isinstance(value, str):
                values[key] = [value]
        return values


class ServerApp(BaseModel):
    verbose: bool
    debug: bool
    host: str
    port: Optional[int]
    open_browser: bool
    force_https: bool
    flask_secret_key: Optional[str]
    generate_cache_control_headers: bool
    server_timing_headers: bool
    csp_directives: Optional[CspDirectives]
    api_base_url: Optional[str]
    web_base_url: Optional[str]

    @root_validator(skip_on_failure=True)
    def check_port(cls, values):
        host = values["host"]
        port = values.get("port")
        if port:
            if not is_port_available(host, port):
                raise ValueError(f"The port selected {port} is in use or invalid, please configure an open port.")
        else:
            values["port"] = find_available_port(host)
        return values

    @root_validator(skip_on_failure=True)
    def check_debug(cls, values):
        if values["debug"]:
            values["verbose"] = True
            values["open_brower"] = False
        else:
            warnings.formatwarning = custom_format_warning
        return values

    @root_validator(skip_on_failure=True)
    def check_api_base_url(cls, values):
        api_base_url = values.get("api_base_url")
        if api_base_url == "local":
            api_base_url = f"http://{values['host']}:{values['port']}"
        if api_base_url and api_base_url.endswith("/"):
            api_base_url = api_base_url[:-1]
        values["api_base_url"] = api_base_url
        return values

    @root_validator(skip_on_failure=True)
    def check_web_base_url(cls, values):
        web_base_url = values["web_base_url"]
        if web_base_url is None:
            web_base_url = values["api_base_url"]
        if web_base_url:
            if web_base_url == "local":
                web_base_url = f"http://{values['host']}:{values['port']}"
            elif web_base_url.endswith("/"):
                web_base_url = web_base_url[:-1]
        values["web_base_url"] = web_base_url
        return values

    @validator("verbose")
    def check_verbose(cls, value):
        if not value:
            sys.tracebacklimit = 0
        else:
            sys.tracebacklimit = 1000
        return value

    @validator("csp_directives")
    def check_csp_directives(cls, value):
        return value if value else {}


class DatarootValue(BaseModel):
    base_url: str
    dataroot: str

    @validator("base_url")
    def check_base_url(cls, base_url):
        bad = False
        # sanity check for well formed base urls
        if os.path.normpath(base_url) != base_url:
            bad = True
        else:
            base_url_parts = base_url.split("/")
            if [quote_plus(part) for part in base_url_parts] != base_url_parts:
                bad = True
            if ".." in base_url_parts:
                bad = True
        if bad:
            raise ValueError(f"invalid base_url={base_url}")
        return base_url


class MultiDataset(BaseModel):
    dataroot: Optional[str] = None
    dataroots: Optional[Dict[str, DatarootValue]] = {}
    index: Union[bool, str] = Field(default=False)

    @root_validator(skip_on_failure=True)
    def check_dataroot(cls, values):
        if all([values["dataroot"], values["dataroots"]]):
            raise ValueError("Must set dataroot or dataroots.")
        elif values["dataroot"]:
            default = dict(base_url="d", dataroot=values["dataroot"])
            values["dataroots"]["d"] = DatarootValue(**default)
            values["dataroot"] = None

        # verify all the base_urls are unique
        base_urls = [d.base_url for d in values["dataroots"].values()]
        if len(base_urls) > len(set(base_urls)):
            raise ValueError("error in multi_dataset__dataroot:  base_urls must be unique")
        # TODO check that at least one dataroot is set. Then we can remove AppConfig.handle_data_source.
        return values


class DataLocator(BaseModel):
    api_base: Optional[str]
    s3_region_name: Union[bool, str, None]


class GeneInfo(BaseModel):
    api_base: Optional[str]


class TiledbCtx(BaseModel):
    sm_tile_cache_size: int = Field(..., alias="sm.tile_cache_size")
    py_init_buffer_bytes: int = Field(..., alias="py.init_buffer_bytes")
    vfs_s3_region: Optional[str] = Field(default=None, alias="vfs.s3.region")


class CxgAdaptor(BaseModel):
    tiledb_ctx: TiledbCtx


class Adaptor(BaseModel):
    cxg_adaptor: CxgAdaptor


class Limits(BaseModel):
    column_request_max: Optional[int]
    diffexp_cellcount_max: Optional[int]


class Server(BaseModel):
    app: ServerApp
    multi_dataset: MultiDataset
    data_locator: DataLocator
    gene_info: GeneInfo
    adaptor: Adaptor
    limits: Limits

    @root_validator(skip_on_failure=True)
    def check_data_locator(cls, values):
        if values["data_locator"].s3_region_name is True:
            path = values["multi_dataset"].dataroots or values["multi_dataset"].dataroot
            # except KeyError as ex:
            #     return values
            if isinstance(path, dict):
                # if multi_dataset.dataroot is a dict, then use the first key
                # that is in s3.   NOTE:  it is not supported to have dataroots
                # in different regions.
                paths = [val.dataroot for val in path.values()]
                for path in paths:
                    if path.startswith("s3://"):
                        break
            if isinstance(path, str) and path.startswith("s3://"):
                region_name = discover_s3_region_name(path)
                if region_name is None:
                    raise ValueError(f"Unable to discover s3 region name from {path}")
                values["data_locator"].s3_region_name = region_name
            else:
                values["data_locator"].s3_region_name = None
        return values

    @root_validator(skip_on_failure=True)
    def check_cxg_adaptor(cls, values):
        if not values["adaptor"].cxg_adaptor.tiledb_ctx.vfs_s3_region and isinstance(
            values["data_locator"].s3_region_name, str
        ):
            values["adaptor"].cxg_adaptor.tiledb_ctx.vfs_s3_region = values["data_locator"].s3_region_name
        return values


class ScriptsItem(BaseModel):
    src: str

    class Config:
        extra = Extra.allow


class DatasetApp(BaseModel):
    scripts: List[Union[str, ScriptsItem]]
    inline_scripts: List[str]
    about_legal_tos: Optional[str]
    about_legal_privacy: Optional[str]

    @validator("scripts")
    def check_scripts(cls, value):
        scripts = []
        for script in value:
            if isinstance(script, str):
                scripts.append(ScriptsItem(src=script))
            else:
                scripts.append(script)
        return scripts


class Presentation(BaseModel):
    max_categories: int
    custom_colors: bool


class Embeddings(BaseModel):
    names: List[str]


class Diffexp(BaseModel):
    enable: bool
    lfc_cutoff: float
    top_n: int
    count: int = 15


class XApproximateDistributionEnum(str, Enum):
    normal = "normal"
    count = "count"


class DefaultDataset(BaseModel):
    app: DatasetApp
    presentation: Presentation
    embeddings: Embeddings
    diffexp: Diffexp
    X_approximate_distribution: XApproximateDistributionEnum = Field(default=XApproximateDistributionEnum.normal)


class AppConfigModel(BaseModel):
    server: Optional[Server]
    default_dataset: Optional[DefaultDataset]
