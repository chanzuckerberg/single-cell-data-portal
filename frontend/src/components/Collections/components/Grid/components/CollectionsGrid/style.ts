import styled from "@emotion/styled";
import Grid from "src/components/common/Grid";
import { VIEW_MODE } from "src/common/hooks/useViewMode";

interface Props {
  mode: VIEW_MODE;
}

const generateGridTemplateColumns = (props: Props) =>
  props.mode === VIEW_MODE.DEFAULT
    ? "minmax(0, 8fr) minmax(0, 6fr) repeat(2, minmax(0, 5fr)) minmax(0, 3fr)"
    : "minmax(0, 65fr) minmax(88px, 24fr) repeat(5, minmax(0, 30fr)) minmax(0, 20fr)";

export const CollectionsGrid = styled(Grid)<Props>`
  grid-template-columns: ${generateGridTemplateColumns};

  /* Collections heading */

  th:first-of-type {
    padding: 0;
  }
`;
