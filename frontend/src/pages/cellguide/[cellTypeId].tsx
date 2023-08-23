import { GetServerSideProps, InferGetServerSidePropsType } from "next";
import CellGuideCard from "src/views/CellGuide/components/CellGuideCard";
import { fetchGptSeoDescription } from "src/common/queries/cellGuide";

const Page = ({
  seoDescription,
  name,
}: InferGetServerSidePropsType<typeof getServerSideProps>): JSX.Element => (
  <CellGuideCard key={name} name={name} seoDescription={seoDescription} />
);

export const getServerSideProps: GetServerSideProps<{
  seoDescription: string;
  name: string;
}> = async (context) => {
  const { params } = context;
  const { cellTypeId: rawCellTypeId } = params ?? {};
  const cellTypeId = (rawCellTypeId as string)?.replace("_", ":");
  const cellType = await fetchGptSeoDescription(cellTypeId);

  const { name, description: seoDescription } = cellType ?? {
    name: "",
    description: "",
  };

  return { props: { seoDescription, name } };
};

export default Page;
