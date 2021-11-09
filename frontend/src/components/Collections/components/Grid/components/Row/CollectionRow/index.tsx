import { Button, Intent, Tag } from "@blueprintjs/core";
import loadable from "@loadable/component";
import Link from "next/link";
import { useRouter } from "next/router";
import { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import { ACCESS_TYPE, VISIBILITY_TYPE } from "src/common/entities";
import {
  useCollection,
  useCreateRevision,
} from "src/common/queries/collections";
import { isTombstonedCollection } from "src/common/utils/typeGuards";
import { aggregateDatasetsMetadata } from "../../../common/utils";
import {
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

const conditionalPopover = (values: string[]) => {
  if (!values || values.length === 0) {
    return <LeftAlignedDetailsCell>-</LeftAlignedDetailsCell>;
  }

  return <AsyncPopover values={values} />;
};

const CollectionRow: FC<Props> = (props) => {
  const { data: collection } = useCollection({
    id: props.id,
    visibility: props.visibility,
  });

  const router = useRouter();

  const navigateToRevision = () => {
    router.push(ROUTES.PRIVATE_COLLECTION.replace(":id", id));
  };

  const [mutate, { isLoading }] = useCreateRevision(navigateToRevision);

  if (!collection || isTombstonedCollection(collection)) return null;

  const handleRevisionClick = () => {
    if (collection?.has_revision === false) {
      mutate(id);
    } else {
      navigateToRevision();
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
        <Link
          href={`/collections/${id}${isPrivate ? "/private" : ""}`}
          passHref
        >
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
            {props.revisionsEnabled && collection.has_revision && (
              <Tag minimal intent={Intent.PRIMARY} data-test-id="revision-tag">
                Revision Pending
              </Tag>
            )}
          </TagContainer>
        )}
      </StyledCell>
      {conditionalPopover(tissue)}
      {conditionalPopover(assay)}
      {conditionalPopover(disease)}
      {conditionalPopover(organism)}
      <RightAlignedDetailsCell>{cell_count || "-"}</RightAlignedDetailsCell>
      {props.revisionsEnabled && visibility === VISIBILITY_TYPE.PUBLIC ? (
        <RevisionCell
          isRevision={collection.has_revision}
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
  isRevision,
  handleRevisionClick,
  isLoading,
}: {
  isRevision?: boolean;
  handleRevisionClick: () => void;
  isLoading: boolean;
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
        {isRevision ? "Continue" : "Start Revision"}
      </Button>
    </RightAlignedDetailsCell>
  );
};

export default CollectionRow;
