import { Position } from "@blueprintjs/core";
import React, { useMemo } from "react";
import RawMoreDropdown from "src/components/common/MoreDropdown";
import { Props as ChooserProps } from "src/components/DropboxChooser/index";
import Menu from "./components/Menu";

interface Props {
  collectionId?: string;
  datasetId?: string;
  revisionsEnabled: boolean;
  onUploadFile: ChooserProps["onUploadFile"];
}

const MoreDropdown = ({
  collectionId = "",
  datasetId = "",
  revisionsEnabled,
  onUploadFile,
}: Props) => {
  const popoverProps = useMemo(() => {
    return {
      content: (
        <Menu
          collectionId={collectionId}
          datasetId={datasetId}
          revisionsEnabled={revisionsEnabled}
          onUploadFile={onUploadFile}
        />
      ),
      position: Position.BOTTOM,
    };
  }, [collectionId, datasetId, revisionsEnabled]);

  return <RawMoreDropdown popoverProps={popoverProps} />;
};

export default MoreDropdown;
