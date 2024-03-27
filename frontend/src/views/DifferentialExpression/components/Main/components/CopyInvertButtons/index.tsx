import { useContext, useMemo } from "react";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import useProcessedQueryGroupFilterDimensions, {
  FilterOptionDimensions,
} from "../common/query_group_filter_dimensions";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import { selectQueryGroup2Filters } from "src/views/DifferentialExpression/common/store/actions";
import { QUERY_GROUP_KEY_TO_FILTER_DIMENSION_MAP } from "../common/constants";

const QUERY_GROUP_KEYS = [
  "tissues",
  "cellTypes",
  "publicationCitations",
  "diseases",
  "ethnicities",
  "sexes",
];
function CopyInvertButtons(): JSX.Element {
  const {
    queryGroups: { queryGroup1, queryGroup2 },
  } = useContext(StateContext);
  const dispatch = useContext(DispatchContext);

  const availableQueryGroup1Filters =
    useProcessedQueryGroupFilterDimensions(queryGroup1);

  const availableQueryGroup2Filters =
    useProcessedQueryGroupFilterDimensions(queryGroup2);

  const handlers = useMemo(() => {
    const handleCopyButtonClickFactory = (key: keyof QueryGroup) => {
      return () => {
        if (!dispatch) {
          return;
        }
        const options = queryGroup1[key];
        const filterDimensionKey = QUERY_GROUP_KEY_TO_FILTER_DIMENSION_MAP[key];
        const optionsWithNames = availableQueryGroup1Filters[
          filterDimensionKey as keyof FilterOptionDimensions
        ].filter((option) => options.includes(option.id));
        dispatch(selectQueryGroup2Filters(key, optionsWithNames));
      };
    };
    const handleInvertButtonClickFactory = (key: keyof QueryGroup) => {
      return () => {
        if (!dispatch) {
          return;
        }
        const options = queryGroup1[key];
        const filterDimensionKey = QUERY_GROUP_KEY_TO_FILTER_DIMENSION_MAP[key];
        const optionsWithNames = availableQueryGroup2Filters[
          filterDimensionKey as keyof FilterOptionDimensions
        ].filter((option) => !options.includes(option.id));

        dispatch(selectQueryGroup2Filters(key, optionsWithNames));
      };
    };

    return QUERY_GROUP_KEYS.map((key) => ({
      key,
      copy: handleCopyButtonClickFactory(key as keyof QueryGroup),
      invert: handleInvertButtonClickFactory(key as keyof QueryGroup),
    }));
  }, [
    dispatch,
    queryGroup1,
    availableQueryGroup1Filters,
    availableQueryGroup2Filters,
  ]);

  return (
    <div>
      {handlers.map(({ key, copy, invert }) => {
        return (
          <div key={key}>
            <button onClick={copy}>Copy {key}</button>
            <button onClick={invert}>Invert {key}</button>
          </div>
        );
      })}
    </div>
  );
}
export default CopyInvertButtons;
