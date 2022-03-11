import Layout from "src/components/Layout";
import { SidebarMainWrapper } from "../style";

interface Props {
  children: JSX.Element;
  SidebarMainWrapperComponent?: typeof SidebarMainWrapper;
}

const SidebarLayout = ({
  children,
  SidebarMainWrapperComponent = SidebarMainWrapper,
}: Props): JSX.Element => {
  return (
    <Layout>
      <SidebarMainWrapperComponent>
        <main>{children}</main>
      </SidebarMainWrapperComponent>
    </Layout>
  );
};

export default SidebarLayout;
