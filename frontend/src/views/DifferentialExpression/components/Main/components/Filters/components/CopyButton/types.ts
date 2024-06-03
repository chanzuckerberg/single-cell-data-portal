import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";

export interface Props {
  queryGroupKey: keyof QueryGroup;
  testId?: string;
}
