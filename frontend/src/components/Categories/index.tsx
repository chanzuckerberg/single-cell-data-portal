// App dependencies
import { FC } from "react";
// Core dependencies
import { CategoryView, OnFilterFn } from "src/common/hooks/useFacetedFilter";
import { Wrapper } from "./style";

interface Props {
  categories: CategoryView[];
  onFilter: OnFilterFn;
}

/**
 * Render categories, category values and their counts. Handle click on category values.
 * @param categories - List of categories to display.
 * @param onFilter - Function to invoke on select of category value.
 */
const Categories: FC<Props> = ({ categories, onFilter }) => {
  return (
    <Wrapper>
      {categories.map((category: CategoryView) => (
        <div key={category.key}>
          <b>{category.key}</b>
          {category.values.map((value) => (
            <div
              key={value.key}
              onClick={() => onFilter(category.key, value.key)}
            >
              {value.selected ? "\u26A1" : null} {value.key} {value.count}
            </div>
          ))}
        </div>
      ))}
    </Wrapper>
  );
};

export default Categories;
