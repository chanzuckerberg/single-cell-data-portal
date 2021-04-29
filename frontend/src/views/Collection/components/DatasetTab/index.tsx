// DEBUG
// DEBUG
// DEBUG
// DEBUG
/* eslint-disable sonarjs/no-duplicate-string */
import { Button, Intent, UL } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import memoize from "lodash/memoize";
import React, { FC, useState } from "react";
import { useQueryCache } from "react-query";
import { Collection, Dataset } from "src/common/entities";
import {
  useCollection,
  useCollectionUploadLinks,
  USE_COLLECTION,
} from "src/common/queries/collections";
import DatasetsGrid from "src/components/Collections/components/Grid/components/DatasetsGrid";
import DropboxChooser, { UploadingFile } from "src/components/DropboxChooser";
import { StyledLink } from "src/views/Collection/common/style";
import { UploadedFiles } from "src/views/Collection/components/ActionButtons";
import DatasetUploadToast from "src/views/Collection/components/DatasetUploadToast";
import EmptyModal from "../EmptyModal";

interface Props {
  collectionID: Collection["id"];
  visibility: Collection["visibility"];
  datasets: Array<Dataset>;
}

const DatasetTab: FC<Props> = ({ collectionID, visibility, datasets }) => {
  // DEBUG
  // DEBUG
  // DEBUG
  // DEBUG
  // DEBUG
  // DEBUG

  datasets = [
    {
      assay: [
        {
          label: "10X 3' v2 sequencing",
          ontology_term_id: "EFO:0009899",
        },
      ],
      cell_count: 5500,
      collection_id: "26ed4278-f88a-4423-ae0a-ad4b90fbcf0d",
      collection_visibility: "PRIVATE",
      created_at: 1619722267.86447,
      dataset_assets: [
        {
          created_at: 1619722362.932458,
          // eslint-disable-next-line @typescript-eslint/ban-ts-comment
          // @ts-ignore
          dataset_id: "65b2656d-9dbc-4f19-bd46-ce720034399e",
          filename: "local.h5ad",
          // eslint-disable-next-line @typescript-eslint/ban-ts-comment
          // @ts-ignore
          filetype: "H5AD",
          id: "c2d0351b-8719-469f-97bc-22a9db599213",
          s3_uri:
            "s3://corpora-data-staging/65b2656d-9dbc-4f19-bd46-ce720034399e/local.h5ad",
          type: "REMIX",
          updated_at: 1619722362.932464,
          user_submitted: true,
        },
        {
          created_at: 1619722371.187526,
          dataset_id: "65b2656d-9dbc-4f19-bd46-ce720034399e",
          filename: "local.loom",
          // eslint-disable-next-line @typescript-eslint/ban-ts-comment
          // @ts-ignore
          filetype: "LOOM",
          id: "fb4f589e-0f05-452a-8199-806866e0f541",
          s3_uri:
            "s3://corpora-data-staging/65b2656d-9dbc-4f19-bd46-ce720034399e/local.loom",
          // eslint-disable-next-line @typescript-eslint/ban-ts-comment
          // @ts-ignore
          type: "REMIX",
          updated_at: 1619722371.187531,
          user_submitted: true,
        },
        {
          created_at: 1619722404.651116,
          dataset_id: "65b2656d-9dbc-4f19-bd46-ce720034399e",
          filename: "local.rds",
          // eslint-disable-next-line @typescript-eslint/ban-ts-comment
          // @ts-ignore
          filetype: "RDS",
          id: "a953e4da-55c6-4506-974c-9967fd6a583e",
          s3_uri:
            "s3://corpora-data-staging/65b2656d-9dbc-4f19-bd46-ce720034399e/local.rds",
          // eslint-disable-next-line @typescript-eslint/ban-ts-comment
          // @ts-ignore
          type: "REMIX",
          updated_at: 1619722404.651123,
          user_submitted: true,
        },
      ],
      dataset_deployments: [
        {
          // eslint-disable-next-line @typescript-eslint/ban-ts-comment
          // @ts-ignore
          created_at: 1619722359.19262,
          dataset_id: "65b2656d-9dbc-4f19-bd46-ce720034399e",
          id: "e1d2e3ba-30bb-416c-a656-2f0f840ba9a5",
          updated_at: 1619722359.192628,
          url:
            "https://cellxgene.staging.single-cell.czi.technology/e/65b2656d-9dbc-4f19-bd46-ce720034399e.cxg/",
        },
      ],
      development_stage: [
        {
          label: "unknown",
          ontology_term_id: "",
        },
      ],
      disease: [
        {
          label: "Alzheimer disease",
          ontology_term_id: "MONDO:0004975",
        },
      ],
      // eslint-disable-next-line @typescript-eslint/ban-ts-comment
      // @ts-ignore
      ethnicity: [
        {
          label: "unknown",
          ontology_term_id: "",
        },
      ],
      id: "65b2656d-9dbc-4f19-bd46-ce720034399e",
      is_valid: false,
      linked_genesets: [],
      name:
        "Molecular characterization of selectively vulnerable neurons in Alzheimer’s Disease: EC astrocytes",
      organism: {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
      processing_status: {
        // eslint-disable-next-line @typescript-eslint/ban-ts-comment
        // @ts-ignore
        conversion_anndata_status: "CONVERTED",
        // eslint-disable-next-line @typescript-eslint/ban-ts-comment
        // @ts-ignore
        conversion_cxg_status: "CONVERTED",
        // eslint-disable-next-line @typescript-eslint/ban-ts-comment
        // @ts-ignore
        conversion_loom_status: "CONVERTED",
        // eslint-disable-next-line @typescript-eslint/ban-ts-comment
        // @ts-ignore
        conversion_rds_status: "CONVERTED",
        created_at: 1619722267.866859,
        dataset_id: "65b2656d-9dbc-4f19-bd46-ce720034399e",
        id: "32cd2d51-0d60-4b3a-a620-15c73b301714",
        // eslint-disable-next-line @typescript-eslint/ban-ts-comment
        // @ts-ignore
        processing_status: "SUCCESS",
        updated_at: 1619722404.706159,
        upload_progress: 1,
        // eslint-disable-next-line @typescript-eslint/ban-ts-comment
        // @ts-ignore
        upload_status: "UPLOADED",

        // eslint-disable-next-line @typescript-eslint/ban-ts-comment
        // @ts-ignore
        validation_status: "VALID",
      },
      published: false,
      revision: 0,
      sex: ["male"],
      tissue: [
        {
          label: "entorhinal cortex",
          ontology_term_id: "UBERON:0002728",
        },
      ],
      updated_at: 1619722361.869618,
    },
  ];

  const CLI_README_LINK =
    "https://github.com/chanzuckerberg/cellxgene/blob/main/dev_docs/schema_guide.md";

  const [uploadLink] = useCollectionUploadLinks(collectionID, visibility);
  const [uploadedFiles, setUploadedFiles] = useState({} as UploadedFiles);
  const { data: collection } = useCollection({ id: collectionID, visibility });

  const queryCache = useQueryCache();

  const isDatasetPresent =
    datasets?.length > 0 || Object.keys(uploadedFiles).length > 0;

  const invalidateCollectionQuery = memoize(
    () => {
      queryCache.invalidateQueries([USE_COLLECTION, collectionID, visibility]);
    },
    () => collectionID + visibility
  );

  const addNewFile = (newFile: UploadingFile) => {
    if (!newFile.link) return;

    const payload = JSON.stringify({ url: newFile.link });

    uploadLink(
      { collectionId: collectionID, payload },
      {
        onSuccess: (datasetID: Dataset["id"]) => {
          newFile.id = datasetID;
          DatasetUploadToast.show({
            icon: IconNames.TICK,
            intent: Intent.PRIMARY,
            message:
              "Your file is being uploaded which will continue in the background, even if you close this window.",
          });
          setUploadedFiles({ ...uploadedFiles, [newFile.id]: newFile });
          queryCache.invalidateQueries(USE_COLLECTION);
        },
      }
    );
  };

  return (
    <>
      {isDatasetPresent ? (
        <DatasetsGrid
          visibility={visibility}
          accessType={collection?.access_type}
          datasets={datasets}
          uploadedFiles={uploadedFiles}
          invalidateCollectionQuery={invalidateCollectionQuery}
        />
      ) : (
        <EmptyModal
          title="No datasets uploaded"
          content={
            <div>
              Before you begin uploading dataset files:
              <UL>
                <li>
                  You must validate your dataset locally. We provide a local CLI
                  script to do this.{" "}
                  <StyledLink href={CLI_README_LINK}>Learn More</StyledLink>
                </li>
                <li>
                  We only support adding datasets in the h5ad format at this
                  time.
                </li>
              </UL>
            </div>
          }
          button={
            <DropboxChooser onUploadFile={addNewFile}>
              <Button
                intent={Intent.PRIMARY}
                outlined
                text={"Add Dataset from Dropbox"}
              />
            </DropboxChooser>
          }
        />
      )}
    </>
  );
};
export default DatasetTab;
