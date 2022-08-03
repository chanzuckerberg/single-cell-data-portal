import Link from "next/link";
import { ROUTES } from "src/common/constants/routes";
import { Details, Header } from "./style";

export default function GetStarted(): JSX.Element {
  return (
    <>
      <Header>Getting Started</Header>
      <Details>
        Use the Add Tissue and Add Gene buttons to find where genes are
        expressed, powered by data from the{" "}
        <Link href={ROUTES.COLLECTIONS} passHref>
          <a href="passHref" rel="noopener" target="_blank">
            cellxgene Data Portal
          </a>
        </Link>
        .
      </Details>
    </>
  );
}
