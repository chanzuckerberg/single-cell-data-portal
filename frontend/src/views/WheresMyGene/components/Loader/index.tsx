import { LoadingIndicator } from "@czi-sds/components";
import { Wrapper } from "./style";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/SaveExport";

export default function Loader(): JSX.Element {
  return (
    <Wrapper
      className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
      data-testid="loading-spinner"
    >
      <LoadingIndicator sdsStyle="tag" />
    </Wrapper>
  );
}
