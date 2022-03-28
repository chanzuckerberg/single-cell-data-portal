import { Header } from "../../common/style";
import { Content, Wrapper } from "./style";

export default function Methodology(): JSX.Element {
  return (
    <Wrapper>
      <Header>Methodology</Header>
      <Content>
        After filtering cells with low coverage (less than 500 genes expressed),
        gene counts are normalized independently for each cell by converting
        counts to quantiles and obtaining the corresponding values from the
        standard normal distribution. Then normalized cell vectors are
        concatenated along the gene axis. The algorithm is summarized in detail
        in our documentation.
      </Content>
    </Wrapper>
  );
}
