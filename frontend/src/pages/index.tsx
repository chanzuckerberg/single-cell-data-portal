import SidebarLayout from "src/components/Layout/components/sidebarLayout";
import Homepage from "src/views/Homepage";

const Page = (): JSX.Element => <Homepage />;

Page.Layout = SidebarLayout;

export default Page;
