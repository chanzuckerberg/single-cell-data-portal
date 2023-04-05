import { useCreateRevision } from "src/common/queries/collections";
import { useRouter } from "next/router";
import { Collection } from "src/common/entities";
import { ROUTES } from "src/common/constants/routes";
import { ActionButton as Button } from "src/views/Collection/components/ActionButtons/style";

interface Props {
  collection: Collection;
}

export default function StartRevisionButton({
  collection,
}: Props): JSX.Element {
  const router = useRouter();
  const navigateToRevision = (id: Collection["id"]) => {
    router.push(ROUTES.COLLECTION.replace(":id", id));
  };
  const { mutate } = useCreateRevision(navigateToRevision);
  return (
    <Button
      onClick={() => mutate(collection.id)}
      sdsStyle="square"
      sdsType="primary"
    >
      Start Revision
    </Button>
  );
}
