import { FEATURES } from "src/common/featureFlags/features";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import DefaultLayout from "src/components/Layout/components/defaultLayout";
import SidebarLayout from "src/components/Layout/components/sidebarLayout";
import Datasets from "src/views/Datasets";

const Page = (): JSX.Element => <Datasets />;

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
