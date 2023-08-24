/* eslint-disable @typescript-eslint/no-empty-interface */
import {
  UseFiltersColumnOptions,
  UseFiltersColumnProps,
  UseFiltersInstanceProps,
  UseFiltersOptions,
  UseFiltersState,
  UseSortByColumnOptions,
  UseSortByColumnProps,
  UseSortByHooks,
  UseSortByInstanceProps,
  UseSortByOptions,
  UseSortByState,
} from "react-table";
import { GridColumnProps } from "src/components/common/Grid/common/entities";

declare module "react-table" {
  export interface TableOptions<D extends Record<string, unknown>>
    extends UseFiltersOptions<D>,
      UseSortByOptions<D> {}

  export interface Hooks<
    D extends Record<string, unknown> = Record<string, unknown>,
  > extends UseSortByHooks<D> {}

  export interface TableInstance<
    D extends Record<string, unknown> = Record<string, unknown>,
  > extends UseFiltersInstanceProps<D>,
      UseSortByInstanceProps<D> {}

  export interface TableState<
    D extends Record<string, unknown> = Record<string, unknown>,
  > extends UseFiltersState<D>,
      UseSortByState<D> {}

  export interface ColumnInterface<
    D extends Record<string, unknown> = Record<string, unknown>,
  > extends UseFiltersColumnOptions<D>,
      UseSortByColumnOptions<D>,
      GridColumnProps {}

  export interface ColumnInstance<
    D extends Record<string, unknown> = Record<string, unknown>,
  > extends UseFiltersColumnProps<D>,
      UseSortByColumnProps<D>,
      GridColumnProps {}
}
