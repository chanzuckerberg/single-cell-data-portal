import styled from "@emotion/styled";
import { CommonThemeProps, fontBodyS, getColors } from "@czi-sds/components";
import { spacesS, spacesXxxs } from "src/common/theme";

interface Props extends CommonThemeProps {
  isOpen: boolean;
}

const grey500 = (props: CommonThemeProps) => getColors(props)?.gray[500];

export const CategoryButton = styled("div")<Props>`
  margin: ${spacesXxxs}px 0 ${spacesS}px;

  .MuiButton-root {
    ${fontBodyS}
    color: ${(props) => (props.isOpen ? `#000000` : grey500(props))};
    font-weight: 500;
    padding: 0;
    text-transform: capitalize;
  }
`;
