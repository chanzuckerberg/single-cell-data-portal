import { EXTERNAL_LINKS } from "src/common/constants/routes";
import { Header } from "../../common/style";
import { Content, Wrapper } from "./style";

export default function Methodology(): JSX.Element {
  return (
    <Wrapper>
      <Header>Methodology</Header>
      <Content>
        After filtering cells with low coverage (less than 500 genes expressed),
        gene counts are normalized independently for each cell by converting
        counts to quantiles and obtaining the corresponding values from a normal
        distribution. Then normalized cell vectors are concatenated along the
        gene axis. The algorithm is summarized in detail in{" "}
        <a href={EXTERNAL_LINKS.WMG_DOC} rel="noopener" target="_blank">
          our documentation
        </a>
        .
      </Content>
    </Wrapper>
  );
}
