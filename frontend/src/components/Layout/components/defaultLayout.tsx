import { useRouter } from "next/router";
import { ReactNode } from "react";
import { ROUTES } from "src/common/constants/routes";
import Layout from "src/components/Layout";
import { DefaultMainWrapper } from "../style";

interface Props {
  children: ReactNode;
}

const DefaultLayout = ({ children }: Props): JSX.Element => {
  const { pathname } = useRouter();

  // REMOVES DEFAULT MAIN WRAPPER ON LANDING PAGE WHICH DISABLES STICKY NAV
  if (pathname === ROUTES.HOMEPAGE) {
    return (
      <Layout>
        <main>{children}</main>
      </Layout>
    );
  } else {
    return (
      <Layout>
        <DefaultMainWrapper
          isCellGuide={pathname.startsWith(ROUTES.CELL_GUIDE)}
        >
          <main>{children}</main>
        </DefaultMainWrapper>
      </Layout>
    );
  }
};

export default DefaultLayout;
