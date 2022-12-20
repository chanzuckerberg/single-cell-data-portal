import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../GeneSearchBar/components/SaveImage";
import Loader from "../Loader";
import { StyledDiv } from "./style";

export default function ScreenTint(): JSX.Element {
  return (
    <>
      <Loader />
      <StyledDiv className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME} />
    </>
  );
}
