import { Position } from "@blueprintjs/core";
import { useMemo } from "react";
import RawMoreDropdown from "src/components/common/MoreDropdown";
import { Props as ChooserProps } from "src/components/DropboxChooser";
import Menu from "./components/Menu";
import { Collection } from "src/common/entities";

interface Props {
  collectionId: Collection["id"];
  datasetId?: string;
  isPublished: boolean;
  revisionsEnabled: boolean;
  onUploadFile: ChooserProps["onUploadFile"];
  isLoading: boolean;
  disabled: boolean;
}

const MoreDropdown = ({
  collectionId,
  datasetId = "",
  isPublished,
  revisionsEnabled,
  onUploadFile,
  isLoading,
  disabled,
}: Props): JSX.Element => {
  const popoverProps = useMemo(() => {
    return {
      content: (
        <Menu
          collectionId={collectionId}
          datasetId={datasetId}
          isPublished={isPublished}
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
    isPublished,
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
