import React from "react";
import { useRouter } from "next/router";
import { CellCardsView } from "./style";

export default function CellCard(): JSX.Element {
  const router = useRouter();

  return (
    <>
      <CellCardsView></CellCardsView>
    </>
  );
}
