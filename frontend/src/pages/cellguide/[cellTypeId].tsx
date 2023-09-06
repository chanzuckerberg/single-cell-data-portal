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
  const cellType = await fetchGptSeoDescription(rawCellTypeId as string);
  const cellTypeMetadata = await fetchCellTypeMetadata();
  const cellTypeId = (rawCellTypeId as string).replace("_", ":");
  const synonyms = cellTypeMetadata[cellTypeId].synonyms;
  const { name, description: seoDescription } = cellType ?? {
    name: "",
    description: "",
  };

  return { props: { seoDescription, name, synonyms } };
};

export default Page;
