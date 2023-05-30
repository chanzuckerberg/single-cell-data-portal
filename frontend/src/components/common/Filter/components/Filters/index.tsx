import React, { Fragment } from "react";
import {
  CategoryView,
  OnFilterFn,
} from "src/components/common/Filter/common/entities";
import Filter from "src/components/common/Filter";
import { FilterDivider } from "./style";

interface Props {
  filters: CategoryView[][];
  onFilter: OnFilterFn;
}

export default function Filters({ filters, onFilter }: Props): JSX.Element {
  return (
    <>
      {filters.map((categoryViews, i) => (
        <Fragment key={i}>
          {i !== 0 && <FilterDivider />}
          <Filter categoryViews={categoryViews} onFilter={onFilter} />
        </Fragment>
      ))}
    </>
  );
}
