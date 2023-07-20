import { GetServerSideProps, InferGetServerSidePropsType } from "next";
import CellGuideCard from "src/views/CellGuide/components/CellGuideCard";
import allCellTypeDescriptionsSEO from "src/views/CellGuide/common/fixtures/allCellTypeDescriptionsSEO.json";

interface AllCellTypeDescriptionsSEO {
  [cellTypeId: string]: { name: string; description: string };
}

const Page = ({
  seoDescription,
  name,
}: InferGetServerSidePropsType<typeof getServerSideProps>): JSX.Element => (
  <CellGuideCard name={name} seoDescription={seoDescription} />
);

export const getServerSideProps: GetServerSideProps<{
  seoDescription: string;
  name: string;
}> = async (context) => {
  const { params } = context;
  const { cellTypeId: rawCellTypeId } = params ?? {};
  const cellTypeId = (rawCellTypeId as string)?.replace("_", ":");

  const cellType = (allCellTypeDescriptionsSEO as AllCellTypeDescriptionsSEO)[
    cellTypeId
  ];

  const { name, description: seoDescription } = cellType;

  return { props: { seoDescription, name } };
};

export default Page;
