import { FC, ReactNode } from "react";
import { DialogLoader } from "src/components/Datasets/components/DownloadDataset/style";
import { FormLabel, FormControl } from "@mui/material";
import { NoneSelected } from "./style";

interface Props {
  downloadPreview?: ReactNode;
  hasDownloadLinks: boolean;
  selected: boolean;
  isLoading: boolean;
}

const Details: FC<Props> = ({
  downloadPreview,
  hasDownloadLinks,
  selected = false,
  isLoading = false,
}) => {
  function renderContent() {
    if (isLoading) {
      return <DialogLoader sdsStyle="minimal" />;
    }

    if (!isLoading && !hasDownloadLinks) {
      return (
        <NoneSelected>
          <h4>No Download Links Available</h4>
          <p>Please try again later.</p>
        </NoneSelected>
      );
    }

    if (!selected) {
      return (
        <NoneSelected>
          <h4>No File Selected</h4>
          <p>Select files to download</p>
        </NoneSelected>
      );
    }

    return downloadPreview;
  }

  return (
    <FormControl>
      <FormLabel>Download Details</FormLabel>
      {renderContent()}
    </FormControl>
  );
};

export default Details;
