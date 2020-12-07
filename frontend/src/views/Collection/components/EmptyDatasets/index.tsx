import { Button, H4, Intent, UL } from "@blueprintjs/core";
import React, { FC } from "react";
import DropboxChooser, {
  Props as DropboxChooserProps,
} from "src/components/DropboxChooser";
import { StyledLink } from "src/views/Collection/common/style";
import { CenterAlignedDiv } from "./style";

const CLI_README_LINK =
  "https://github.com/chanzuckerberg/cellxgene/blob/main/dev_docs/schema_guide.md";

interface Props {
  onSelectUploadLink: DropboxChooserProps["onSelectUploadLink"];
}

const EmptyDatasets: FC<Props> = ({ onSelectUploadLink }) => {
  return (
    <CenterAlignedDiv>
      <H4>No datasets uploaded</H4>
      <div>
        Before you begin uploading dataset files:
        <UL>
          <li>
            You must validate your dataset locally. We provide a local CLI
            script to do this.{" "}
            <StyledLink href={CLI_README_LINK}>Learn More</StyledLink>
          </li>
          <li>
            We only support adding datasets in the h5ad format at this time.
          </li>
        </UL>
      </div>
      <DropboxChooser onSelectUploadLink={onSelectUploadLink}>
        <Button
          intent={Intent.PRIMARY}
          outlined
          text={"Add Dataset from Dropbox"}
        />
      </DropboxChooser>
    </CenterAlignedDiv>
  );
};

export default EmptyDatasets;
