/**
 * IMPORTANT(thuang): `frontend/src/pages/cellguide/tissues/[tissueId]/cell-types/[cellTypeId].tsx`
 * imports the exports from this file, since both routes use the same components.
 * If we export more things from this file, we need to export them from the other file as well,
 * to make sure that /cellguide/tissues/[tissueId]/cell-types/[cellTypeId] route continue
 * to work.
 */

import { GetServerSideProps, InferGetServerSidePropsType } from "next";
import CellGuideCard from "src/views/CellGuide/components/CellGuideCard";
import {
  fetchGptSeoDescription,
  fetchCellTypeMetadata,
} from "src/common/queries/cellGuide";

const Page = ({
  seoDescription,
  name,
  synonyms,
}: InferGetServerSidePropsType<typeof getServerSideProps>): JSX.Element => (
  <CellGuideCard
    key={name}
    name={name}
    seoDescription={seoDescription}
    synonyms={synonyms}
  />
);

export const getServerSideProps: GetServerSideProps<{
  seoDescription: string;
  name: string;
  synonyms?: string[];
}> = async (context) => {
  const { params } = context;
  const { cellTypeId: rawCellTypeId } = params ?? {};
  const seoDescription = await fetchGptSeoDescription(rawCellTypeId as string);
  const cellTypeMetadata = await fetchCellTypeMetadata();
  const cellTypeId = (rawCellTypeId as string).replace("_", ":");
  const { synonyms, name } = cellTypeMetadata[cellTypeId];
  return { props: { seoDescription, name, synonyms } };
};

export default Page;
