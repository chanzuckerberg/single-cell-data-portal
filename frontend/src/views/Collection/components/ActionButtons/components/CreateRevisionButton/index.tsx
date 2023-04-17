import { ActionButton as Button } from "src/views/Collection/components/ActionButtons/style";

interface Props {
  handleCreateRevision: () => void;
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
