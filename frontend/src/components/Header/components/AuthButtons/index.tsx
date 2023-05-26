import { AnchorButton, Button, MenuDivider } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { Menu, MenuItem } from "@czi-sds/components";
import { MouseEventHandler, useState } from "react";
import { API } from "src/common/API";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useUserInfo } from "src/common/queries/auth";
import { AuthButtonWrapper } from "src/components/Header/style";
import { API_URL } from "src/configs/configs";
import CuratorAPIKeyGenerator from "../CuratorAPIKeyGenerator";
import { LogOutAnchor, LogOutEmail, LogOutText } from "./style";

const AuthButtons = (): JSX.Element | null => {
  const hasAuth = get(FEATURES.CURATOR) === BOOLEAN.TRUE;

  const { data: userInfo, isLoading, error } = useUserInfo(hasAuth);

  if (userInfo && error) {
    // (thuang): Force refresh page to log user out
    window.location.reload();
  }

  if (!hasAuth || isLoading) return null;

  return (
    <AuthButtonWrapper>
      {userInfo?.name ? (
        <LoggedInButtons name={userInfo.name} email={userInfo.email} />
      ) : (
        <LoggedOutButtons />
      )}
    </AuthButtonWrapper>
  );
};

function LoggedInButtons({ name, email }: { name?: string; email?: string }) {
  const authName = isEmail(name) ? "Account" : name;
  const [anchorEl, setAnchorEl] = useState<Element | null>(null);

  const handleClick: MouseEventHandler<HTMLElement> = (event) => {
    setAnchorEl(event.currentTarget);
  };

  const handleClose = () => {
    setAnchorEl(null);
  };

  return (
    <>
      <Content />
      <Button
        onClick={handleClick}
        minimal
        rightIcon={IconNames.CHEVRON_DOWN}
        text={authName}
      />
    </>
  );

  function Content() {
    const curatorAPIFeature = get(FEATURES.CURATOR);
    return (
      <Menu
        anchorEl={anchorEl}
        keepMounted
        open={Boolean(anchorEl)}
        onClose={handleClose}
      >
        {curatorAPIFeature && (
          <div>
            <CuratorAPIKeyGenerator />
            <MenuDivider />
          </div>
        )}
        <LogOutAnchor href={`${API_URL}${API.LOG_OUT}`}>
          <MenuItem data-testid="log-out">
            <LogOutText>Logout</LogOutText>
            <LogOutEmail data-testid="user-email">{email}</LogOutEmail>
          </MenuItem>
        </LogOutAnchor>
      </Menu>
    );
  }
}

function LoggedOutButtons() {
  return (
    <AnchorButton
      href={`${API_URL}${API.LOG_IN}`}
      minimal
      rightIcon={IconNames.CHEVRON_DOWN}
      text="Log In"
    />
  );
}

function isEmail(name?: string): boolean {
  if (!name) return false;

  return name.includes("@");
}

export default AuthButtons;
