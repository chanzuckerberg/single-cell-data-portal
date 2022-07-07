import { Button, Intent, Tag } from "@blueprintjs/core";
import loadable from "@loadable/component";
import Link from "next/link";
import { useRouter } from "next/router";
import { FC } from "react";
import { PLURALIZED_METADATA_LABEL } from "src/common/constants/metadata";
import { ROUTES } from "src/common/constants/routes";
import { ACCESS_TYPE, Collection, VISIBILITY_TYPE } from "src/common/entities";
import {
  useCollection,
  useCreateRevision,
} from "src/common/queries/collections";
import { isTombstonedCollection } from "src/common/utils/typeGuards";
import { aggregateDatasetsMetadata } from "../../../common/utils";
import {
  DiseaseDetailsCell,
  LeftAlignedDetailsCell,
  RightAlignedDetailsCell,
  StyledCell,
  StyledRow,
  TagContainer,
} from "../common/style";
import { CollectionTitleText } from "./style";

interface Props {
  id: string;
  accessType?: ACCESS_TYPE;
  visibility: VISIBILITY_TYPE;
  revisionsEnabled?: boolean;
}

const AsyncPopover = loadable(
  () =>
    /*webpackChunkName: 'CollectionRow/components/Popover' */ import(
      "src/components/Collections/components/Grid/components/Popover"
    )
);

const AsyncDiseasePopover = loadable(
  () =>
    /*webpackChunkName: 'src/components/common/Grid/components/DiseaseCell' */ import(
      "src/components/common/Grid/components/DiseaseCell"
    )
);

const conditionalPopover = (
  label: PLURALIZED_METADATA_LABEL,
  values: string[]
) => {
  if (!values || values.length === 0) {
    return <LeftAlignedDetailsCell>-</LeftAlignedDetailsCell>;
  }

  return <AsyncPopover label={label} values={values} />;
};

const conditionalDiseasePopover = (
  label: PLURALIZED_METADATA_LABEL,
  values: string[]
) => {
  if (!values || values.length === 0) {
    return <LeftAlignedDetailsCell>-</LeftAlignedDetailsCell>;
  }

  return (
    <DiseaseDetailsCell>
      <AsyncDiseasePopover label={label} values={values} />
    </DiseaseDetailsCell>
  );
};

const CollectionRow: FC<Props> = (props) => {
  const { data: collection } = useCollection({
    id: props.id,
  });

  const router = useRouter();

  const navigateToRevision = (id: Collection["id"]) => {
    router.push(ROUTES.COLLECTION.replace(":id", id));
  };

  const { mutate, isLoading } = useCreateRevision(navigateToRevision);

  if (!collection || isTombstonedCollection(collection)) return null;

  const handleRevisionClick = () => {
    if (!collection?.revisioning_in) {
      mutate(id);
    } else {
      navigateToRevision(collection?.revisioning_in || collection.id);
    }
  };

  // If there is an explicity accessType only show collections with that accessType
  if (props.accessType && collection.access_type !== props.accessType) {
    return null;
  }

  const { id, visibility, name } = collection;
  const datasets = Array.from(collection.datasets.values());
  const isPrivate = visibility === VISIBILITY_TYPE.PRIVATE;

  const { tissue, assay, disease, organism, cell_count } =
    aggregateDatasetsMetadata(datasets);

  return (
    <StyledRow data-test-id="collection-row">
      <StyledCell>
        <Link href={`/collections/${id}`} passHref>
          <CollectionTitleText data-test-id="collection-link" href="passHref">
            {name}
          </CollectionTitleText>
        </Link>

        {props.accessType === ACCESS_TYPE.WRITE && (
          <TagContainer>
            <Tag
              minimal
              intent={isPrivate ? Intent.PRIMARY : Intent.SUCCESS}
              data-test-id="visibility-tag"
            >
              {isPrivate ? "Private" : "Published"}
            </Tag>
            {props.revisionsEnabled && collection.revisioning_in && (
              <Tag minimal intent={Intent.PRIMARY} data-test-id="revision-tag">
                Revision Pending
              </Tag>
            )}
          </TagContainer>
        )}
      </StyledCell>
      {conditionalPopover(PLURALIZED_METADATA_LABEL.TISSUE, tissue)}
      {conditionalPopover(PLURALIZED_METADATA_LABEL.ASSAY, assay)}
      {conditionalDiseasePopover(PLURALIZED_METADATA_LABEL.DISEASE, disease)}
      {conditionalPopover(PLURALIZED_METADATA_LABEL.ORGANISM, organism)}
      <RightAlignedDetailsCell>{cell_count || "-"}</RightAlignedDetailsCell>
      {props.revisionsEnabled && visibility === VISIBILITY_TYPE.PUBLIC ? (
        <RevisionCell
          revisionId={collection.revisioning_in}
          handleRevisionClick={handleRevisionClick}
          isLoading={isLoading}
        />
      ) : (
        <RightAlignedDetailsCell />
      )}
    </StyledRow>
  );
};

const RevisionCell = ({
  handleRevisionClick,
  isLoading,
  revisionId,
}: {
  handleRevisionClick: () => void;
  isLoading: boolean;
  revisionId: Collection["revisioning_in"];
}) => {
  return (
    <RightAlignedDetailsCell>
      <Button
        loading={isLoading}
        intent={Intent.PRIMARY}
        minimal
        onClick={handleRevisionClick}
        data-test-id="revision-action-button"
      >
        {revisionId ? "Continue" : "Start Revision"}
      </Button>
    </RightAlignedDetailsCell>
  );
};

export default CollectionRow;
