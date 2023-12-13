export interface Props {
  id?: string;
  isHubSpotReady?: boolean;
  setError: (error: string) => void;
  setEmail: (email: string) => void;
  email: string;
  emailValidationError?: string;
}
