/**
 * Helpers for schema management
 *
 * TODO: all this would be much more natural if done with a framework
 * like immutable.js
 */

import fromEntries from "../fromEntries";
import { RawSchema, Schema } from "../../common/types/schema";

/**
 * System wide schema assumptions:
 *  - schema and data wil be consistent (eg, for user-created annotations)
 *  - schema will be internally self-consistent (eg, index matches columns)
 */

/**
 * Index schema for ease of use
 */
export function indexEntireSchema(schema: RawSchema): Schema {
  (schema as Schema).annotations.obsByName = fromEntries(
    schema.annotations?.obs?.columns?.map((v) => [v.name, v]) || [],
  );
  (schema as Schema).annotations.varByName = fromEntries(
    schema.annotations?.var?.columns?.map((v) => [v.name, v]) || [],
  );
  (schema as Schema).layout.obsByName = fromEntries(
    schema.layout?.obs?.map((v) => [v.name, v]) || [],
  );
  (schema as Schema).layout.varByName = fromEntries(
    schema.layout?.var?.map((v) => [v.name, v]) || [],
  );

  return schema as Schema;
}
