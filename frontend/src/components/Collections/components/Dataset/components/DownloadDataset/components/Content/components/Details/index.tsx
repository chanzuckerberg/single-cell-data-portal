import { FC, ReactNode } from "react";
import { DialogLoader } from "src/components/Datasets/components/DownloadDataset/style";
import { FormLabel } from "@mui/material";
import {
  FormControl,
  SeuratNotice,
  StyledIcon,
  StyledLink,
  TextWrapper,
} from "./style";

export const PROMPT_TEXT =
  "Select one of the data formats to view its download details.";

interface Props {
  downloadPreview?: ReactNode;
  selected: boolean;
  fileSize: number;
  isLoading: boolean;
  selectedFormat: string;
}

const MEGA_BYTES = 2 ** 20;

const DOC_SITE_URL = "/docs/03__Download%20Published%20Data#seurat-deprecated";

const Details: FC<Props> = ({
  downloadPreview,
  selected = false,
  fileSize = 0,
  isLoading = false,
  selectedFormat,
}) => {
  function renderContent() {
    if (isLoading) {
      return <DialogLoader sdsStyle="minimal" />;
    }

    if (!selected) {
      return <div>{PROMPT_TEXT}</div>;
    }

    return <div>{`${Math.round(fileSize / MEGA_BYTES)}MB`}</div>;
  }

  return (
    <FormControl>
      {/* TODO: Remove after Seurat has been deprecated */}
      {selectedFormat == "RDS" && (
        <SeuratNotice>
          <StyledIcon
            sdsIcon="ExclamationMarkCircle"
            sdsSize="l"
            sdsType="static"
          />
          <TextWrapper>
            Seurat support will be removed between Nov 15 - Dec 31, 2024. You
            can download and convert the .h5ad yourself by following these {""}
            <StyledLink href={DOC_SITE_URL}>instructions</StyledLink>.
          </TextWrapper>
        </SeuratNotice>
      )}
      <FormLabel>Download Details</FormLabel>
      {renderContent()}
      {downloadPreview}
    </FormControl>
  );
};

export default Details;
