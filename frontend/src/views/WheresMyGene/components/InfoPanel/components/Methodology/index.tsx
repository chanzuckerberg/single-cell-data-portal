import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ROUTES } from "src/common/constants/routes";
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
        <a
          href={`${ROUTES.DOCS}/04__Analyze%20Public%20Data/4_2__Gene%20Expression%20Documentation/4_2_1__Get%20Started`}
          rel="noopener"
          target="_blank"
          onClick={() => {
            track(EVENTS.WMG_SOURCE_DOCUMENTATION_CLICKED);
          }}
        >
          our documentation
        </a>
        .
      </Content>
    </Wrapper>
  );
}
