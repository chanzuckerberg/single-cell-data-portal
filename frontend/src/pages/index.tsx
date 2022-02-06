import { FEATURES } from "src/common/featureFlags/features";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import DefaultLayout from "src/components/Layout/components/defaultLayout";
import SidebarLayout from "src/components/Layout/components/sidebarLayout";
import Homepage from "src/views/Homepage";

const Page = (): JSX.Element => <Homepage />;

Page.Layout = function _Layout({
  children,
}: {
  children: JSX.Element;
}): JSX.Element {
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  const Layout = isFilterEnabled ? SidebarLayout : DefaultLayout;

  return <Layout>{children}</Layout>;
};

export default Page;
