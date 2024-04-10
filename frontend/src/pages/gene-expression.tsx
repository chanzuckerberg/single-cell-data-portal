import SidebarLayout from "src/components/Layout/components/sidebarLayout";
import { SideBarLayoutMainWrapper } from "src/views/WheresMyGeneV2/style";
import WheresMyGeneV2 from "src/views/WheresMyGeneV2";

const Page = (): JSX.Element => <WheresMyGeneV2 />;

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
