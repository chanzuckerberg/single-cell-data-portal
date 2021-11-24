/* eslint-disable @typescript-eslint/no-empty-interface */
import {
  UseFiltersColumnOptions,
  UseFiltersColumnProps,
  UseFiltersInstanceProps,
  UseFiltersOptions,
  UseFiltersState,
} from "react-table";

declare module "react-table" {
  export interface TableOptions<D extends Record<string, unknown>>
    extends UseFiltersOptions<D> {}

  export interface TableInstance<
    D extends Record<string, unknown> = Record<string, unknown>
  > extends UseFiltersInstanceProps<D> {}

  export interface TableState<
    D extends Record<string, unknown> = Record<string, unknown>
  > extends UseFiltersState<D> {}

  export interface ColumnInterface<
    D extends Record<string, unknown> = Record<string, unknown>
  > extends UseFiltersColumnOptions<D> {}

  export interface ColumnInstance<
    D extends Record<string, unknown> = Record<string, unknown>
  > extends UseFiltersColumnProps<D> {}
}
