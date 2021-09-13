import { Position } from "@blueprintjs/core";
import { useMemo } from "react";
import RawMoreDropdown from "src/components/common/MoreDropdown";
import { Props as ChooserProps } from "src/components/DropboxChooser/index";
import Menu from "./components/Menu";

interface Props {
  collectionId?: string;
  datasetId?: string;
  revisionsEnabled: boolean;
  onUploadFile: ChooserProps["onUploadFile"];
  isLoading: boolean;
  disabled: boolean;
}

const MoreDropdown = ({
  collectionId = "",
  datasetId = "",
  revisionsEnabled,
  onUploadFile,
  isLoading,
  disabled,
}: Props) => {
  const popoverProps = useMemo(() => {
    return {
      content: (
        <Menu
          collectionId={collectionId}
          datasetId={datasetId}
          revisionsEnabled={revisionsEnabled}
          onUploadFile={onUploadFile}
          isLoading={isLoading}
        />
      ),
      disabled,
      position: Position.BOTTOM,
    };
  }, [
    collectionId,
    datasetId,
    revisionsEnabled,
    isLoading,
    onUploadFile,
    disabled,
  ]);

  return (
    <RawMoreDropdown popoverProps={popoverProps} buttonProps={{ disabled }} />
  );
};

export default MoreDropdown;
