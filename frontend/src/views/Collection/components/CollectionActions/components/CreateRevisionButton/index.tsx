import { ActionButton as Button } from "src/views/Collection/components/CollectionActions/style";
import { CreateRevisionFn } from "src/views/Collection/components/CollectionActions";

interface Props {
  handleCreateRevision: CreateRevisionFn;
}

export default function CreateRevisionButton({
  handleCreateRevision,
}: Props): JSX.Element {
  return (
    <Button onClick={handleCreateRevision} sdsStyle="square" sdsType="primary">
      Start Revision
    </Button>
  );
}
