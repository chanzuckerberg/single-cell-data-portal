import React, {  } from "react";
import { useRouter } from "next/router";
import {
  CellCardsView,
} from "./style";

// enum of available descriptions
type DescriptionOptions = "GPT3.5" | "Wikipedia" | "OLS v4";
const availableDescriptions: DescriptionOptions[] = [
  "GPT3.5",
  "Wikipedia",
  "OLS v4",
];

export default function CellCard(): JSX.Element {
  const router = useRouter();

  return (
    <>
      <CellCardsView>

      </CellCardsView>
    </>
  );
}
