import {
  AnchorButton,
  Button,
  Menu,
  MenuItem,
  Popover,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { API } from "src/common/API";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useUserInfo } from "src/common/queries/auth";
import { API_URL } from "src/configs/configs";

const AuthButtons = (): JSX.Element | null => {
  const hasAuth = get(FEATURES.CURATOR) === BOOLEAN.TRUE;

  const { data: userInfo, isLoading, error } = useUserInfo(hasAuth);

  if (userInfo && error) {
    // (thuang): Force refresh page to log user out
    window.location.reload();
  }

  if (!hasAuth || isLoading) return null;

  return userInfo?.name ? (
    <LoggedInButtons name={userInfo.name} email={userInfo.email} />
  ) : (
    <LoggedOutButtons />
  );
};

function LoggedInButtons({ name, email }: { name?: string; email?: string }) {
  const authName = isEmail(name) ? "Account" : name;

  return (
    <Popover content={<Content />}>
      <Button minimal rightIcon={IconNames.CHEVRON_DOWN} text={authName} />
    </Popover>
  );

  function Content() {
    return (
      <Menu>
        <MenuItem data-testid="user-email" text={`Logged in as: ${email}`} />
        <MenuItem
          data-testid="log-out"
          text="Log Out"
          href={`${API_URL}${API.LOG_OUT}`}
          icon={IconNames.LOG_OUT}
        />
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
