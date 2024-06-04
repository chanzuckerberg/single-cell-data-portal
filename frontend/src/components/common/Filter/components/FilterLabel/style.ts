import styled from "@emotion/styled";
import { CommonThemeProps, fontBodyS, getColors } from "@czi-sds/components";
import { spacesS, spacesXxs, spacesXxxs } from "src/common/theme";

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
    gap: ${spacesXxs}px;
    justify-content: flex-start;
    letter-spacing: -0.006em;
    min-width: 0;
    padding: 0;
    text-transform: capitalize;

    &:hover {
      color: #000000;
    }

    .MuiButton-endIcon {
      align-items: center;
      display: flex;
      height: 16px;
      justify-content: center;
      margin: 0;
      width: 16px;
    }
  }
`;
