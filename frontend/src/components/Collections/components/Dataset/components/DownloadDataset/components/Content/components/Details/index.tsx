import { FC, ReactNode } from "react";
import { DialogLoader } from "src/components/Datasets/components/DownloadDataset/style";
import { FormLabel } from "@mui/material";
import { FormControl } from "./style";

export const PROMPT_TEXT =
  "Select one of the data formats to view its download details.";

interface Props {
  curlPreview?: ReactNode;
  selected: boolean;
  fileSize: number;
  isLoading: boolean;
}

const MEGA_BYTES = 2 ** 20;

const Details: FC<Props> = ({
  curlPreview,
  selected = false,
  fileSize = 0,
  isLoading = false,
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
      <FormLabel>Download Details</FormLabel>
      {renderContent()}
      {curlPreview}
    </FormControl>
  );
};

export default Details;
