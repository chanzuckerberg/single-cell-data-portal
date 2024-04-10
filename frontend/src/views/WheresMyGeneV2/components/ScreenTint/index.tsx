import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/SaveExport";
import Loader from "../Loader";
import { StyledDiv } from "./style";

export default function ScreenTint({
  isDownloading,
}: {
  isDownloading: {
    isLoading: boolean;
    blur?: boolean;
  };
}): JSX.Element {
  return (
    <>
      {isDownloading.isLoading && (
        <>
          <Loader />
          <StyledDiv
            className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
            style={{
              backdropFilter: isDownloading.blur ? "blur(10px)" : "none",
            }}
          />
        </>
      )}
    </>
  );
}
