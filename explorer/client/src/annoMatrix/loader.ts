import { doBinaryRequest, doFetch } from "./fetchHelpers";
import { matrixFBSToDataframe } from "../util/stateManager/matrix";
import { _getColumnSchema } from "./schema";
import { _whereCacheCreate, WhereCache } from "./whereCache";
import AnnoMatrix from "./annoMatrix";
import PromiseLimit from "../util/promiseLimit";
import {
  _expectComplexQuery,
  _expectSimpleQuery,
  _hashStringValues,
  _urlEncodeComplexQuery,
  _urlEncodeLabelQuery,
  _urlOptionalEncodeNbinsSuffix,
  ComplexQuery,
  Query,
} from "./query";
import { normalizeResponse } from "./normalize";
import { Field, RawSchema } from "../common/types/schema";
import { Dataframe } from "../util/dataframe";

const promiseThrottle = new PromiseLimit<ArrayBuffer>(5);

export default class AnnoMatrixLoader extends AnnoMatrix {
  baseURL: string;

  /*
  AnnoMatrix implementation which proxies to HTTP server using the CXG REST API.
  Used as the base (non-view) instance.

  Public API is same as AnnoMatrix class (refer there for API description),
  with the addition of the constructor which bootstraps:

    new AnnoMatrixLoader(serverBaseURL, schema) -> instance

  */
  constructor(baseURL: string, schema: RawSchema) {
    const { nObs, nVar } = schema.dataframe;
    super(schema, nObs, nVar);

    if (baseURL[baseURL.length - 1] !== "/") {
      // must have trailing slash
      baseURL += "/";
    }
    this.baseURL = baseURL;
    Object.seal(this);
  }

  /**
   ** Private below
   **/
  async _doLoad(
    field: Field,
    query: Query,
    nBins: number | null,
  ): Promise<[WhereCache | null, Dataframe]> {
    /*
    _doLoad - evaluates the query against the field. Returns:
      * whereCache update: column query map mapping the query to the column labels
      * Dataframe containing the new columns (one per dimension)
    */
    let doRequest;
    let priority = 10; // default fetch priority

    switch (field) {
      case "obs":
      case "var": {
        doRequest = _obsOrVarLoader(this.baseURL, field, query, nBins);
        break;
      }
      case "X": {
        doRequest = _XLoader(this.baseURL, field, query, nBins);
        break;
      }
      case "emb": {
        doRequest = _embLoader(this.baseURL, field, query, nBins);
        priority = 0; // high prio load for embeddings
        break;
      }
      default:
        throw new Error("Unknown field name");
    }
    const buffer = await promiseThrottle.priorityAdd(priority, doRequest);
    let result = matrixFBSToDataframe(buffer);
    if (!result || result.isEmpty()) throw Error("Unknown field/col");

    const whereCacheUpdate = _whereCacheCreate(
      field,
      query,
      result.colIndex.labels(),
    );

    result = normalizeResponse(field, this.schema, result);

    return [whereCacheUpdate, result];
  }
}

/*
Utility functions below
*/

function _embLoader(
  baseURL: string,
  _field: Field,
  query: Query,
  nBins: number | null = null,
): () => Promise<ArrayBuffer> {
  _expectSimpleQuery(query);

  const urlBase = `${baseURL}layout/obs`;
  const urlQuery = _urlEncodeLabelQuery("layout-name", query);
  const url = _urlOptionalEncodeNbinsSuffix(`${urlBase}?${urlQuery}`, nBins);
  return () => doBinaryRequest(url);
}

function _obsOrVarLoader(
  baseURL: string,
  field: Field,
  query: Query,
  nBins: number | null = null,
): () => Promise<ArrayBuffer> {
  _expectSimpleQuery(query);

  const urlBase = `${baseURL}annotations/${field}`;
  const urlQuery = _urlEncodeLabelQuery("annotation-name", query);
  const url = _urlOptionalEncodeNbinsSuffix(`${urlBase}?${urlQuery}`, nBins);
  return () => doBinaryRequest(url);
}

function _XLoader(
  baseURL: string,
  _field: Field,
  query: Query,
  nBins: number | null = null,
): () => Promise<ArrayBuffer> {
  _expectComplexQuery(query);

  // Casting here as query is validated to be complex in _expectComplexQuery above.
  const complexQuery = query as ComplexQuery;

  if ("where" in complexQuery) {
    const urlBase = `${baseURL}data/var`;
    const urlQuery = _urlEncodeComplexQuery(complexQuery);
    const url = _urlOptionalEncodeNbinsSuffix(`${urlBase}?${urlQuery}`, nBins);
    return () => doBinaryRequest(url);
  }

  if ("summarize" in complexQuery) {
    const urlBase = `${baseURL}summarize/var`;
    const urlQuery = _urlEncodeComplexQuery(complexQuery);

    if (urlBase.length + urlQuery.length < 2000) {
      const url = _urlOptionalEncodeNbinsSuffix(
        `${urlBase}?${urlQuery}`,
        nBins,
      );
      return () => doBinaryRequest(url);
    }

    const url = _urlOptionalEncodeNbinsSuffix(
      `${urlBase}?key=${_hashStringValues([urlQuery])}`,
      nBins,
    );
    return async () => {
      const res = await doFetch(url, {
        method: "POST",
        body: urlQuery,
        headers: new Headers({
          Accept: "application/octet-stream",
          "Content-Type": "application/x-www-form-urlencoded",
        }),
      });
      return res.arrayBuffer();
    };
  }

  throw new Error("Unknown query structure");
}
