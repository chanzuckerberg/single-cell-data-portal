import SidebarLayout from "src/components/Layout/components/sidebarLayout";
import WheresMyGene from "src/views/WheresMyGene";
import { SideBarLayoutMainWrapper } from "src/views/WheresMyGene/style";

const Page = (): JSX.Element => <WheresMyGene />;

// This layout needs to be defined in order to correctly style the SideBars and their wrappers
Page.Layout = function Layout({
  children,
}: {
  children: JSX.Element;
}): JSX.Element {
  return (
    <SidebarLayout SidebarMainWrapperComponent={SideBarLayoutMainWrapper}>
      {children}
    </SidebarLayout>
  );
};

export default Page;
