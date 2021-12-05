import SidebarLayout from "src/components/Layout/components/sidebarLayout";
import Datasets from "src/views/Datasets";

const Page = (): JSX.Element => <Datasets />;

Page.Layout = SidebarLayout;

export default Page;
