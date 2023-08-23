import { GetServerSideProps, InferGetServerSidePropsType } from "next";
import { fetchTissueMetadata } from "src/common/queries/cellGuide";
import TissueCard from "src/views/CellGuide/components/TissueCard";

const Page = ({
  description,
  name,
}: InferGetServerSidePropsType<typeof getServerSideProps>): JSX.Element => (
  <TissueCard key={name} name={name} description={description} />
);

export const getServerSideProps: GetServerSideProps<{
  description: string;
  name: string;
}> = async (context) => {
  const { params } = context;
  const { tissueId: rawTissueId } = params ?? {};
  const allTissues = await fetchTissueMetadata();
  const tissueEntry = allTissues[rawTissueId as string];
  const name = tissueEntry?.name || "";
  const description = tissueEntry?.uberonDescription || "";

  return { props: { description, name } };
};

export default Page;
