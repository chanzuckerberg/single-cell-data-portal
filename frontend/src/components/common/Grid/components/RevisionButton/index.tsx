import React from "react";
import { ROUTES } from "src/common/constants/routes";
import { useRouter } from "next/router";
import { RevisionButton as IconButton } from "./style";

interface Props {
  id?: string;
}

export default function RevisionButton({ id }: Props): JSX.Element | null {
  const router = useRouter();
  return id ? (
    <IconButton
      data-testid="view-collection-revision"
      sdsIcon="open"
      sdsSize="small"
      sdsType="primary"
      onClick={() => router.push(ROUTES.COLLECTION.replace(":id", id))}
    />
  ) : null;
}
