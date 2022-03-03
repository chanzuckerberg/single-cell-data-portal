import SidebarLayout from "src/components/Layout/components/sidebarLayout";
import WheresMyGene from "src/views/WheresMyGene";

const Page = (): JSX.Element => <WheresMyGene />;

// This layout needs to be defined in order to correctly style the SideBars and their wrappers
Page.Layout = SidebarLayout;

export default Page;
