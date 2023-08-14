import { GetServerSideProps, InferGetServerSidePropsType } from "next";
import TissueCard from "src/views/CellGuide/components/TissueCard";
import { allTissues, allTissueDescriptions } from "src/pages/api/tissue_cards";

interface Tissues {
  id: string;
  label: string;
}

interface AllTissueDescriptions {
  [tissueId: string]: string;
}

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
  const tissueId = (rawTissueId as string)?.replace("_", ":");

  const name =
    (allTissues as Tissues[]).find(
      (tissue: { id: string; label: string }) => tissue.id === tissueId
    )?.label || "";

  const description = (allTissueDescriptions as AllTissueDescriptions)[
    tissueId
  ];

  return { props: { description, name } };
};

export default Page;
