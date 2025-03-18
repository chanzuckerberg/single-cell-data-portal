import { FC, ReactNode } from "react";
import { DialogLoader } from "src/components/Datasets/components/DownloadDataset/style";
import { FormLabel } from "@mui/material";
import {
  FormControl,
  NoneSelected
} from "./style";

interface Props {
  downloadPreview?: ReactNode;
  selected: boolean;
  isLoading: boolean;
}

const Details: FC<Props> = ({
  downloadPreview,
  selected = false,
  isLoading = false,
}) => {
  function renderContent() {
    if (isLoading) {
      return <DialogLoader sdsStyle="minimal" />;
    }

    if (!selected) {
      return <NoneSelected>
        <h4>No File Selected</h4>
        <p>Select files to download</p>
      </NoneSelected>
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
