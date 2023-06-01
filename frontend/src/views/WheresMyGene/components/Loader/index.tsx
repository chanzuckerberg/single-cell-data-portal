import { LoadingIndicator } from "@czi-sds/components";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../GeneSearchBar/components/SaveExport";
import { Wrapper } from "./style";

export default function Loader(): JSX.Element {
  return (
    <Wrapper className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
      <LoadingIndicator sdsStyle="tag" />
    </Wrapper>
  );
}
